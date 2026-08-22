#!/usr/bin/env python3
"""Test-only strict-ALA thermodynamic closure oracle.

The executable checks reduced algebra and static production contracts.  It does
not run CitcomS finite-element assembly, SUPG, or production MPI communicators;
those evidence gates remain explicitly unimplemented.
"""

from __future__ import annotations

import argparse
from collections import Counter
import configparser
import hashlib
import json
import math
import re
import subprocess
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Iterable

import numpy as np


GLOBAL_ROOT = Path(__file__).resolve().parents[2]
PROJECT_ROOT = Path(__file__).resolve().parents[4]
RUNS_ROOT = PROJECT_ROOT / "runs"
DEFAULT_CONFIG = RUNS_ROOT / "cmbhf_ALA_strict.cfg"
DEFAULT_REFSTATE = RUNS_ROOT / "refstate_ALA_strict.txt"
DEFAULT_INTERVAL = RUNS_ROOT / "interval_ALA_strict.txt"
DEFAULT_MESH = RUNS_ROOT / "GLB.coor.global.dat"
DEFAULT_KATSURA = RUNS_ROOT / "prepare_data" / "Katsura2022_TableS5_thermoelastic.csv"

PHASE_NAMES = ("410", "520", "660")
BOUNDARIES_KM = np.asarray((410.0, 520.0, 660.0))
SERIAL_BUOYANCY_TOLERANCE = 1.0e-12
MPI_BUOYANCY_TOLERANCE = 1.0e-10
TRAJECTORY_INTEGRATION_TOLERANCE_K = 1.0e-7
TRAJECTORY_CLOSURE_TOLERANCE_K = 1.0e-3

TOPOLOGY_COMPONENTS = {
    "full": GLOBAL_ROOT / "CitcomS" / "Components" / "Sphere" / "FullSphere.py",
    "regional": GLOBAL_ROOT
    / "CitcomS"
    / "Components"
    / "Sphere"
    / "RegionalSphere.py",
}

BETA_SOURCE_ROLES = {
    "beta_ala": {
        "lib/global_defs.h": "reference-state field declaration",
        "lib/Material_properties.c": (
            "strict nodal input candidate plus inactive legacy initialization/diagnostics"
        ),
    },
    "ala_beta_supplied": {
        "lib/global_defs.h": "reference-state field declaration",
        "lib/Material_properties.c": "endpoint-average candidate construction",
        "lib/Stokes_flow_Incomp.c": "explicit causal diagnostic alternative",
    },
    "ala_beta_density": {
        "lib/global_defs.h": "reference-state field declaration",
        "lib/Material_properties.c": "density-log-secant candidate construction",
        "lib/Stokes_flow_Incomp.c": "explicit causal diagnostic alternative",
    },
    "ala_beta_interval": {
        "lib/global_defs.h": "reference-state field declaration",
        "lib/Material_properties.c": "file-backed interval candidate and selection",
    },
    "ala_beta": {
        "lib/global_defs.h": "reference-state field declaration",
        "lib/Material_properties.c": "single selected runtime field",
        "lib/Element_calculations.c": "continuity, transpose, and AL consumers",
        "lib/Drive_solvers.c": "pressure-work diagnostic consumer",
    },
    "gamma_eff": {
        "lib/global_defs.h": "reference-state field declaration",
        "lib/Material_properties.c": "thermodynamic closure diagnostic",
    },
}

# Frozen production-code occurrence multiset. These hashes are SHA-256 over
# whitespace-normalized, comment-masked source lines; they are never derived
# from the source being audited. Counter comparison detects both new uses and
# duplicated/deleted existing uses while allowing unrelated line-number drift.
BETA_SOURCE_OCCURRENCE_CONTRACT = (
    ("ala_beta", "lib/Drive_solvers.c", "5bd7656316e859a7210129278c0f02e2b22ea3f38194020e66e9a85c836b2a17"),
    ("ala_beta", "lib/Element_calculations.c", "9c09aeb7301685da0b50cbc5e2cdd872cfeee20999fa8ecbf1d67311fb6627f3"),
    ("ala_beta", "lib/Element_calculations.c", "9cb7075ceb63446dab74fb1b80a3323a0b939991659e908b1426a3c2102e5d76"),
    ("ala_beta", "lib/Element_calculations.c", "fdbd7c7d1560fe7355687e8bbbc8dca810d463ee15601322d5460aee92125a51"),
    ("ala_beta", "lib/Material_properties.c", "6b67aa37e582741df6e4abd5091f8526f5e3336cd01a1c2baf30df343c261488"),
    ("ala_beta_supplied", "lib/Material_properties.c", "13e14866b9c5e30b1b3770124582d9ea2bc45f32f4961cb5c4b60a160c974a1c"),
    ("ala_beta_density", "lib/Material_properties.c", "da211503bc61c7490e32254f14cc956584ad205cf62b7369059e0cfdd98793f1"),
    ("ala_beta_interval", "lib/Material_properties.c", "f8b95ae343f14a1021deeeebb58f76c655b7e4f25b35e215ef171cb6ba0f0854"),
    ("beta_ala", "lib/Material_properties.c", "8630c7d29284fb57daa76bfa0be1f5dfadd37391536e398e82246332f4cf5f1c"),
    ("ala_beta", "lib/Material_properties.c", "e97de0a126133cfc77e80396eebb077d2fc07fcb3eaa4f8295d75f4545158034"),
    ("ala_beta_supplied", "lib/Material_properties.c", "e56ca8f61a5fabe2f244e85ec38e559921fd285734c3c6c1d51299f4e1a13afa"),
    ("ala_beta_density", "lib/Material_properties.c", "5b540cd9f2a33714a8805c17c0277ee6456fa23f5fd6f72743696b75dfb885f7"),
    ("ala_beta_interval", "lib/Material_properties.c", "e1a325fd137b6183f2d547208065f9e66850eb62e40df04f234c446868e34fcb"),
    ("beta_ala", "lib/Material_properties.c", "3dc2f8969c6366bb0b69ba86482dbb67baa262ea6d475bf42d3269ffd9bff39e"),
    ("gamma_eff", "lib/Material_properties.c", "6d158d90060e2ee0c54c460442a6938056915f2f439e2d2ae2e44af0fdd93ee1"),
    ("gamma_eff", "lib/Material_properties.c", "8da820a002dc2558686dad24f14827f7cd7e6d11b1f68db4a34cd4c309016086"),
    ("ala_beta", "lib/Material_properties.c", "0839e71205098280e4ff9935eac0424d706fcdc090fbf0f2d47d60a61edce93f"),
    ("ala_beta", "lib/Material_properties.c", "1bcdc175fe35416032539fb46b763c85b1a1c2c9a194cb5a3b0f34a18ca1e214"),
    ("ala_beta_supplied", "lib/Material_properties.c", "c477e520788a66eab81ea031caa84a7778440fe5a9210e1dae1d80c00869aefc"),
    ("ala_beta_supplied", "lib/Material_properties.c", "265bc8acec0827f1a568965bd68833b95967f64c3ffccac5016d51ebbaf8ace3"),
    ("ala_beta_density", "lib/Material_properties.c", "6f19b4ff8b39f99b78bab0b2fb1864dbe44bc90e4cabd6595511ff6361febb66"),
    ("ala_beta_density", "lib/Material_properties.c", "85f54b9de6897eb024d03b9b2fa2b6767bc1f85c73669049fb6c4fb0b911cb0f"),
    ("ala_beta_interval", "lib/Material_properties.c", "484cf2a97531d6a7f243b2f15720e0e302a4161a2ea3fa0a8863c40daefd5432"),
    ("ala_beta_interval", "lib/Material_properties.c", "fc897d2dec78fa72ff407f58f4d1d674ef2d4ba95ec527368b45d770a2947aaf"),
    ("beta_ala", "lib/Material_properties.c", "6e966a8c46bae465ab490447d60eb7b0bacfc7d6a46740a017a3e31face45992"),
    ("beta_ala", "lib/Material_properties.c", "c01b94ea1152e038ccb97171666e372e830ffc43e27659fa76ff1ab5b4db31ad"),
    ("beta_ala", "lib/Material_properties.c", "7cb14c4b8d968dcb7b8bfbaa76514bc7c6c21daf160868c4e25b8976524315bc"),
    ("gamma_eff", "lib/Material_properties.c", "7cb14c4b8d968dcb7b8bfbaa76514bc7c6c21daf160868c4e25b8976524315bc"),
    ("ala_beta_supplied", "lib/Material_properties.c", "f128b9f268ab6d5a78c8ea0556101196ae83661d933849f7acb0194ab3dbe108"),
    ("beta_ala", "lib/Material_properties.c", "3352a8c16b7747a9d2cd21869f27a1128f3013c9c50954720766ea01952586a1"),
    ("beta_ala", "lib/Material_properties.c", "6279b56a329ee5b1cca1c90f66e44b9fabe66ec8b49ab1948b3ee57e9f11de93"),
    ("ala_beta_density", "lib/Material_properties.c", "880e59f63f6071cce21e6b80ba8ee9df9cdee7ebb5e74905ca8435c9f7ae0724"),
    ("ala_beta_supplied", "lib/Material_properties.c", "07131079afb7251728695ecb6c7849d07a813017b3d573cd969cee4036e5bf55"),
    ("ala_beta_supplied", "lib/Material_properties.c", "d18ba072d04c028d78d5d9d6b08d3b69c65d06b8b7a16b385c986b76b9648989"),
    ("ala_beta_interval", "lib/Material_properties.c", "168abd1027e8e153945b3e22dda93327d226bb567719379586545e4520a5e618"),
    ("ala_beta_density", "lib/Material_properties.c", "8391450db0e3b21416e7c8fb4bd7e295d6db004a8bd2e9b7f0d80d3d36311109"),
    ("ala_beta_supplied", "lib/Material_properties.c", "0688e71c9a34abbc920430403389d50f43e6922f9a4133f63e9d5b605d4f9371"),
    ("ala_beta_supplied", "lib/Material_properties.c", "ffc4b30d2e96c20ce3c449da0e8614e64f9e650085b6db6450b5196701dea6a3"),
    ("ala_beta_density", "lib/Material_properties.c", "cc924c36f6165b1cf44499e1a189e12c687b0b9cf40b20a0b6f5a5eb786e6b7b"),
    ("ala_beta_density", "lib/Material_properties.c", "815c37766b25a7409eab90511de06320c18a1556e09982a59392a94ca1d82700"),
    ("ala_beta_supplied", "lib/Material_properties.c", "191e285a6649bd69cbf845b64e55ca032800069068320f389482cdd75c5ec138"),
    ("ala_beta_density", "lib/Material_properties.c", "e22bcbe05f5684da9555993dca59391e7601e06a363a617870b5fc64d8efa989"),
    ("ala_beta", "lib/Material_properties.c", "0b7ccb25a4c9c00a6b7d4412f8f0ae8fbb055fc1818f0bdcbcbb68864d79eead"),
    ("gamma_eff", "lib/Material_properties.c", "242b5842ddc25b9967d8dfe1a415fae2f97d6dc4e1f0e439953d253a26b2d556"),
    ("gamma_eff", "lib/Material_properties.c", "a7de59b3c32bbd3fa4bc8ff8884ea4f90a0e016e6c82a5738b3232f4029b21e0"),
    ("ala_beta", "lib/Material_properties.c", "4a7158a16e801fa425ee0937a5a0357e8cdbe0e1ce5b805761ab486deaa1b30d"),
    ("beta_ala", "lib/Material_properties.c", "e07f0b67cd772446a8c7e243fcef2f17532e8bfb9654fd299877b79f940b85c2"),
    ("gamma_eff", "lib/Material_properties.c", "4d8aade9c45496fb445ec3de10f6a484afa8c3c6ea8cd30ca0aa597c31145800"),
    ("beta_ala", "lib/Material_properties.c", "0acf2230530f51f41b8c227607a988b8ab95c36fe81621236bb42db633ef3fe4"),
    ("gamma_eff", "lib/Material_properties.c", "dfc9c48c8489f585e4b7c3b73ad2d0205b04f110e394df015d552d9681af46a5"),
    ("beta_ala", "lib/Material_properties.c", "98b6ebc03ebe0e8f851180e89eceb0a6c0f9dd7a459f36086dea7ee75b685e19"),
    ("gamma_eff", "lib/Material_properties.c", "df78ade5f14a9fd788072037c0f100b0989e07329c00fa5ec9ee1cf43b783b40"),
    ("gamma_eff", "lib/Material_properties.c", "e5fe1c099e9cb05220ddb5bc5910843fe9b5a8fc9eaa24f668aadfbe501a49bc"),
    ("beta_ala", "lib/Material_properties.c", "ceeb85b96b29211f4d60264d612c68702feea4bccb42c008b2b1e7be8e2a7ba0"),
    ("ala_beta_interval", "lib/Material_properties.c", "9ea603ea50418324f53bb5cec88d2b143d5f5d00bddd5509abf265cf7b8f4161"),
    ("ala_beta_interval", "lib/Material_properties.c", "ee18273ff3244a9a2532a9bd688befd81764e5a7ad1e1c0eef56423c5aa63d1d"),
    ("ala_beta_interval", "lib/Material_properties.c", "867e9db81cdadaa9d886319b4eee2274d783cd83c2a567046f981e33396197ce"),
    ("gamma_eff", "lib/Material_properties.c", "019cb59924adaf20444c8d30983dde41c863ac63a0e5efb6c113ae345b2268be"),
    ("beta_ala", "lib/Material_properties.c", "12ab050ee64424b5a8cb808b7394522e8d20a6ced52ea039f6e2ab199b6448a7"),
    ("ala_beta_density", "lib/Stokes_flow_Incomp.c", "b1af900ec7843fe265837229084def328157e296ef804a9d3f774cb51e14f4c5"),
    ("ala_beta_supplied", "lib/Stokes_flow_Incomp.c", "f1d650595066a36455d8dd92768c905264e753d5b11bd26d2352672dc9769a46"),
)

BETA_REF_STATE_DECLARATION_CONTRACT = (
    "has_beta_ala",
    "has_beta_interval",
    "beta_interval_filename",
    "ala_beta",
    "ala_beta_supplied",
    "ala_beta_density",
    "ala_beta_interval",
    "beta_ala",
    "gamma_eff",
)

IMPLEMENTED = "IMPLEMENTED"
NOT_IMPLEMENTED = "NOT-IMPLEMENTED"
RUN = "RUN"
NOT_RUN = "NOT-RUN"
BLOCKED = "BLOCKED"
PASS = "PASS"
FAIL = "FAIL"


@dataclass(frozen=True)
class OracleInputs:
    config: Path = DEFAULT_CONFIG
    refstate: Path = DEFAULT_REFSTATE
    interval: Path = DEFAULT_INTERVAL
    mesh: Path = DEFAULT_MESH
    katsura: Path = DEFAULT_KATSURA


@dataclass(frozen=True)
class EnergyTerms:
    R_storage: np.ndarray
    R_advection: np.ndarray
    R_adiabatic: np.ndarray
    R_phase: np.ndarray
    R_diffusion: np.ndarray
    R_viscous: np.ndarray
    R_internal: np.ndarray
    R_boundary: np.ndarray
    R_stabilization: np.ndarray

    @property
    def total(self) -> np.ndarray:
        values = asdict(self).values()
        return sum((value for value in values), np.zeros_like(self.R_storage))


def _parser(path: Path) -> configparser.ConfigParser:
    parser = configparser.ConfigParser(interpolation=None, strict=True)
    parser.optionxform = str
    if not parser.read(path, encoding="utf-8"):
        raise FileNotFoundError(path)
    return parser


def _vector(section: configparser.SectionProxy, name: str) -> np.ndarray:
    values = np.asarray(
        [float(item.strip()) for item in section[name].split(",")], dtype=float
    )
    if values.shape != (3,) or not np.all(np.isfinite(values)):
        raise ValueError(f"{name} must contain three finite values")
    return values


def _norms(value: np.ndarray) -> dict[str, float]:
    array = np.asarray(value, dtype=float)
    if array.size == 0 or not np.all(np.isfinite(array)):
        raise ValueError("norm input must be nonempty and finite")
    return {
        "L2": float(np.sqrt(np.sum(array * array))),
        "RMS": float(np.sqrt(np.mean(array * array))),
        "Linf": float(np.max(np.abs(array))),
    }


def _relative_norms(residual: dict[str, float], scale: dict[str, float]) -> dict[str, float]:
    return {
        name: residual[name] / max(scale[name], 1.0e-300)
        for name in ("L2", "RMS", "Linf")
    }


def _relative_rms(value: np.ndarray, scale: np.ndarray) -> float:
    numerator = float(np.sqrt(np.mean(np.asarray(value, dtype=float) ** 2)))
    denominator = max(
        float(np.sqrt(np.mean(np.asarray(scale, dtype=float) ** 2))), 1.0e-300
    )
    return numerator / denominator


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _git_metadata(path: Path) -> dict[str, Any]:
    def run(*arguments: str) -> str:
        result = subprocess.run(
            ("git", "-C", str(path), *arguments),
            check=True,
            capture_output=True,
            text=True,
        )
        return result.stdout.strip()

    tracked_diff = subprocess.run(
        ("git", "-C", str(path), "diff", "--binary", "HEAD", "--", "."),
        check=True,
        capture_output=True,
    ).stdout
    return {
        "commit": run("rev-parse", "HEAD"),
        "branch": run("branch", "--show-current"),
        "dirty": bool(run("status", "--porcelain")),
        "tracked_diff_sha256": hashlib.sha256(tracked_diff).hexdigest(),
    }


def _configured_surface_processors(
    cfg: configparser.ConfigParser,
) -> tuple[str, int, str]:
    solver = cfg["CitcomS"].get("solver", "").strip().lower()
    if solver not in TOPOLOGY_COMPONENTS:
        raise ValueError(f"unsupported CitcomS solver topology: {solver!r}")

    mesher = cfg["CitcomS.solver.mesher"]
    explicit = mesher.get("nproc_surf", fallback=None)
    if explicit is not None:
        count = int(explicit)
        provenance = "explicit CitcomS.solver.mesher.nproc_surf"
    else:
        component = TOPOLOGY_COMPONENTS[solver]
        source = component.read_text(encoding="utf-8")
        match = re.search(
            r'nproc_surf\s*=\s*pyre\.inventory\.int\('
            r'\s*["\']nproc_surf["\']\s*,\s*default\s*=\s*(\d+)\s*\)',
            source,
        )
        if match is None:
            raise ValueError(f"cannot resolve nproc_surf default from {component}")
        count = int(match.group(1))
        provenance = f"default from {component.relative_to(GLOBAL_ROOT)}"
    if count < 1:
        raise ValueError("nproc_surf must be positive")
    return solver, count, provenance


def _compiled_library_sources() -> list[Path]:
    makefile = GLOBAL_ROOT / "lib" / "Makefile.am"
    text = makefile.read_text(encoding="utf-8")
    match = re.search(
        r"^sources\s*=\s*\\\n(?P<body>.*?)(?:\n\n|\Z)",
        text,
        re.MULTILINE | re.DOTALL,
    )
    if match is None:
        raise ValueError(f"cannot resolve compiled source inventory from {makefile}")
    names = re.findall(r"\b[A-Za-z0-9_]+\.(?:c|h)\b", match.group("body"))
    paths = [GLOBAL_ROOT / "lib" / name for name in names]
    missing = [path for path in paths if not path.is_file()]
    if missing:
        raise FileNotFoundError(missing[0])
    return paths


_C_NONCODE = re.compile(
    r'//[^\n]*|/\*.*?\*/|"(?:\\.|[^"\\])*"|\'(?:\\.|[^\'\\])*\'',
    re.DOTALL,
)
_BETA_FIELD = re.compile(
    r"(?:E->)?refstate\."
    r"(beta_ala|ala_beta_supplied|ala_beta_density|"
    r"ala_beta_interval|ala_beta|gamma_eff)\b"
)


def _mask_c_noncode(source: str) -> str:
    """Mask C comments and literals while retaining line-number provenance."""

    return _C_NONCODE.sub(
        lambda match: "".join(
            "\n" if character == "\n" else " " for character in match.group()
        ),
        source,
    )


def _normalized_line_hash(line: str) -> str:
    normalized = " ".join(line.split())
    return hashlib.sha256(normalized.encode("utf-8")).hexdigest()


def _counter_records(counter: Counter[tuple[str, str, str]]) -> list[dict[str, Any]]:
    return [
        {
            "field": field,
            "file": relative,
            "normalized_line_sha256": line_hash,
            "count": count,
        }
        for (field, relative, line_hash), count in sorted(counter.items())
    ]


def scan_beta_source_occurrences(
    sources: dict[str, str],
    expected_contract: Iterable[tuple[str, str, str]] = BETA_SOURCE_OCCURRENCE_CONTRACT,
) -> dict[str, Any]:
    """Compare known refstate beta/Gamma accesses with a frozen exact multiset."""

    expected = Counter(expected_contract)
    actual: Counter[tuple[str, str, str]] = Counter()
    occurrences: list[dict[str, Any]] = []
    for relative, source in sources.items():
        masked = _mask_c_noncode(source)
        for line_number, line in enumerate(masked.splitlines(), start=1):
            line_hash = _normalized_line_hash(line)
            for match in _BETA_FIELD.finditer(line):
                field = match.group(1)
                key = (field, relative, line_hash)
                actual[key] += 1
                occurrences.append(
                    {
                        "field": field,
                        "file": relative,
                        "line": line_number,
                        "normalized_line_sha256": line_hash,
                        "allowlisted_role": BETA_SOURCE_ROLES.get(field, {}).get(
                            relative
                        ),
                    }
                )

    unexpected = actual - expected
    missing = expected - actual
    return {
        "expected_occurrence_count": sum(expected.values()),
        "actual_occurrence_count": sum(actual.values()),
        "occurrences": occurrences,
        "unexpected_occurrences": _counter_records(unexpected),
        "missing_occurrences": _counter_records(missing),
        "contract_matches": not unexpected and not missing,
        "allowlist_scope": (
            "frozen field/file/whitespace-normalized-line-hash multiset; "
            "C comments and literals excluded"
        ),
    }


def ref_state_beta_gamma_declarations(source: str) -> tuple[str, ...]:
    """Inventory beta/Gamma-like members declared by struct REF_STATE."""

    masked = _mask_c_noncode(source)
    match = re.search(r"struct\s+REF_STATE\s*\{(?P<body>.*?)\}\s*;", masked, re.DOTALL)
    if match is None:
        raise ValueError("cannot locate struct REF_STATE")
    identifiers = re.findall(r"\b[A-Za-z_][A-Za-z0-9_]*\b", match.group("body"))
    return tuple(
        name for name in identifiers if "beta" in name.lower() or "gamma" in name.lower()
    )


def _phase_fraction(
    depth: np.ndarray,
    rho: np.ndarray,
    gravity: np.ndarray,
    temperature: np.ndarray,
    phase_depth: float,
    clapeyron: float,
    trans_temperature: float,
    inverse_width: float,
) -> np.ndarray:
    pressure = (
        (depth - phase_depth) * rho * gravity
        - clapeyron * (temperature - trans_temperature)
    )
    return 0.5 * (1.0 + np.tanh(inverse_width * pressure))


def horizontal_projection(field: np.ndarray, partitions: int = 1) -> np.ndarray:
    """Unweighted reduced projection; not production FE area weighting."""
    values = np.asarray(field, dtype=float)
    if values.ndim != 2:
        raise ValueError("horizontal projection expects (columns, radii)")
    if partitions < 1 or partitions > values.shape[0]:
        raise ValueError("invalid horizontal partition count")
    sums = np.zeros(values.shape[1], dtype=float)
    counts = np.zeros(values.shape[1], dtype=np.int64)
    for indices in np.array_split(np.arange(values.shape[0]), partitions):
        sums += np.sum(values[indices], axis=0, dtype=float)
        counts += len(indices)
    return values - sums / counts


class StrictALAClosureOracle:
    def __init__(self, inputs: OracleInputs = OracleInputs()) -> None:
        self.inputs = inputs
        self.cfg = _parser(inputs.config)
        self.refstate = np.loadtxt(inputs.refstate, comments="#", ndmin=2)
        self.intervals = np.loadtxt(inputs.interval, comments="#", ndmin=2)
        self.radius = np.loadtxt(inputs.mesh, skiprows=1, usecols=1)
        if self.refstate.shape[1] != 8:
            raise ValueError("phase-inclusive strict refstate must have eight columns")
        if len(self.radius) != len(self.refstate):
            raise ValueError("mesh and refstate row counts differ")
        if self.intervals.shape != (len(self.radius) - 1, 4):
            raise ValueError("interval beta shape differs from radial elements")

        const = self.cfg["CitcomS.solver.const"]
        solver = self.cfg["CitcomS.solver"]
        mesher = self.cfg["CitcomS.solver.mesher"]
        tracer = self.cfg["CitcomS.solver.tracer"]
        self.radius_m = const.getfloat("radius")
        self.delta_temperature_k = const.getfloat("Tbottom") - const.getfloat("Ttop")
        self.surface_temperature = const.getfloat("Ttop") / self.delta_temperature_k
        self.rho0 = const.getfloat("rho0")
        self.cp0 = const.getfloat("Cp0")
        self.g0 = const.getfloat("g0")
        self.alpha0 = const.getfloat("alpha0")
        self.rayleigh = solver.getfloat("rayleigh")
        self.dissipation_number = solver.getfloat("dissipation_number")
        self.nprocx = mesher.getint("nprocx")
        self.nprocy = mesher.getint("nprocy")
        self.nprocz = mesher.getint("nprocz")
        self.topology_solver, self.mpi_caps, self.topology_provenance = (
            _configured_surface_processors(self.cfg)
        )
        self.mpi_ranks = self.mpi_caps * self.nprocx * self.nprocy * self.nprocz
        self.horizontal_comm_size = self.mpi_caps * self.nprocx * self.nprocy
        self._phase_parameters_loaded = False
        self.buoyancy_ratios = np.asarray(
            [float(item.strip()) for item in tracer["buoyancy_ratio"].split(",")],
            dtype=float,
        )
        if not np.all(np.isfinite(self.buoyancy_ratios)):
            raise ValueError("buoyancy_ratio must be finite")

        self.depth = 1.0 - self.radius
        self.rho = self.refstate[:, 0]
        self.gravity = self.refstate[:, 1]
        self.tref = self.refstate[:, 2]
        self.alpha = self.refstate[:, 3]
        self.cp = self.refstate[:, 4]
        if not np.all(np.diff(self.radius) > 0.0):
            raise ValueError("mesh must be ordered CMB to surface")

    def _load_phase_parameters(self) -> None:
        if self._phase_parameters_loaded:
            return
        phase = self.cfg["CitcomS.solver.phase"]
        self.phase_depth = _vector(phase, "phase_depth")
        self.phase_density_jump = _vector(phase, "phase_delta_rho")
        density_scale = self.rho0 * self.alpha0 * self.delta_temperature_k
        self.phase_ra = self.rayleigh * self.phase_density_jump / density_scale
        self.phase_width = _vector(phase, "phase_width")
        self.phase_clapeyron = _vector(phase, "phase_clapeyron")
        self.phase_trans_temperature = _vector(phase, "phase_transT")
        self._phase_parameters_loaded = True

    def input_manifest(self) -> dict[str, Any]:
        files = {
            "configuration": self.inputs.config,
            "reference_state": self.inputs.refstate,
            "interval_beta": self.inputs.interval,
            "mesh": self.inputs.mesh,
            "Katsura_source": self.inputs.katsura,
        }
        inspected_sources = {
            str(path.relative_to(GLOBAL_ROOT)): _sha256(path)
            for path in (
                GLOBAL_ROOT / "lib" / "Advection_diffusion.c",
                GLOBAL_ROOT / "lib" / "Phase_change.c",
                GLOBAL_ROOT / "lib" / "Pan_problem_misc_functions.c",
                GLOBAL_ROOT / "lib" / "Material_properties.c",
                GLOBAL_ROOT / "lib" / "Element_calculations.c",
                GLOBAL_ROOT / "lib" / "Drive_solvers.c",
                GLOBAL_ROOT / "lib" / "Stokes_flow_Incomp.c",
                GLOBAL_ROOT / "lib" / "global_defs.h",
                GLOBAL_ROOT / "lib" / "Makefile.am",
                TOPOLOGY_COMPONENTS[self.topology_solver],
            )
        }
        return {
            "source": _git_metadata(GLOBAL_ROOT),
            "runs": _git_metadata(RUNS_ROOT),
            "files": {
                name: {"path": str(path), "sha256": _sha256(path)}
                for name, path in files.items()
            },
            "test_artifacts": {
                path.name: {"path": str(path), "sha256": _sha256(path)}
                for path in (
                    Path(__file__).resolve(),
                    Path(__file__).with_name(
                        "test_thermodynamic_closure_oracle.py"
                    ),
                    Path(__file__).with_name("mpi_thermodynamic_closure.c"),
                    Path(__file__).with_name(
                        "STAGE0_THERMODYNAMIC_CLOSURE_REPORT.md"
                    ),
                    Path(__file__).with_name("README.md"),
                )
            },
            "inspected_source_sha256": inspected_sources,
            "topology": {
                "caps": self.mpi_caps,
                "nproc_surf": self.mpi_caps,
                "nprocx": self.nprocx,
                "nprocy": self.nprocy,
                "nprocz": self.nprocz,
                "world_ranks": self.mpi_ranks,
                "horizontal_comm_size": self.horizontal_comm_size,
                "horizontal_communicators": self.nprocz,
                "solver": self.topology_solver,
                "provenance": self.topology_provenance,
            },
            "oracle_thresholds": {
                "trajectory_closure_K": TRAJECTORY_CLOSURE_TOLERANCE_K,
                "RK4_refinement_K": TRAJECTORY_INTEGRATION_TOLERANCE_K,
                "provenance": (
                    "ORACLE-DEFINED; the Stage 0 request supplied no numeric "
                    "trajectory closure or integration threshold"
                ),
            },
        }

    def phase_fraction(self, phase: int, temperature: np.ndarray) -> np.ndarray:
        self._load_phase_parameters()
        return _phase_fraction(
            self.depth,
            self.rho,
            self.gravity,
            np.asarray(temperature),
            self.phase_depth[phase],
            self.phase_clapeyron[phase],
            self.phase_trans_temperature[phase],
            1.0 / self.phase_width[phase],
        )

    def temperature_diagnostics(
        self,
        temperature: np.ndarray | None = None,
        *,
        include_phases: bool = True,
    ) -> dict[str, Any]:
        total = self.tref.copy() if temperature is None else np.asarray(temperature)
        phases = None
        if include_phases:
            phases = []
            for index, name in enumerate(PHASE_NAMES):
                phases.append(
                    {
                        "phase": name,
                        "X_i": self.phase_fraction(index, total).tolist(),
                        "X_i_ref": self.phase_fraction(index, self.tref).tolist(),
                    }
                )
        return {
            "T_total": total.tolist(),
            "T_ref": self.tref.tolist(),
            "deltaT": (total - self.tref).tolist(),
            "phases": phases,
            "phase_diagnostics_status": RUN if include_phases else BLOCKED,
        }

    def manufactured_radial_velocity(self, amplitude: float = 1.0) -> np.ndarray:
        return amplitude / (self.radius * self.radius * self.rho)

    def mass_flux_constancy(self, amplitude: float = 1.0) -> dict[str, float]:
        velocity = self.manufactured_radial_velocity(amplitude)
        flux = self.radius * self.radius * self.rho * velocity
        residual = flux - amplitude
        return {
            "maximum_absolute_residual": float(np.max(np.abs(residual))),
            "relative_RMS": _relative_rms(residual, flux),
        }

    def _interpolate(self, field: np.ndarray, radius: float) -> float:
        return float(np.interp(radius, self.radius, field))

    def _material_temperature_gradient(
        self, radius: float, temperature: float, phases: Iterable[int] = ()
    ) -> float:
        rho = self._interpolate(self.rho, radius)
        gravity = self._interpolate(self.gravity, radius)
        alpha = self._interpolate(self.alpha, radius)
        cp = self._interpolate(self.cp, radius)
        theta = temperature + self.surface_temperature
        phase_linear = 0.0
        phase_quadratic = 0.0
        phase_indices = tuple(phases)
        if phase_indices:
            self._load_phase_parameters()
        for index in phase_indices:
            depth = 1.0 - radius
            fraction = float(
                _phase_fraction(
                    np.asarray([depth]),
                    np.asarray([rho]),
                    np.asarray([gravity]),
                    np.asarray([temperature]),
                    self.phase_depth[index],
                    self.phase_clapeyron[index],
                    self.phase_trans_temperature[index],
                    1.0 / self.phase_width[index],
                )[0]
            )
            localization = 2.0 * fraction * (1.0 - fraction) / self.phase_width[index]
            coefficient = self.phase_ra[index] / self.rayleigh
            phase_linear += localization * self.phase_clapeyron[index] * coefficient
            phase_quadratic += (
                localization * self.phase_clapeyron[index] ** 2 * coefficient
            )
        inertia = 1.0 + self.dissipation_number * theta * phase_quadratic
        return -(
            self.dissipation_number
            * alpha
            * gravity
            * theta
            * (1.0 + phase_linear)
            / (cp * inertia)
        )

    def _integrate_segment(
        self, start: float, stop: float, temperature: float, substeps: int
    ) -> tuple[float, list[tuple[float, float]]]:
        step = (stop - start) / substeps
        radius = start
        path = []
        for _ in range(substeps):
            k1 = self._material_temperature_gradient(radius, temperature)
            k2 = self._material_temperature_gradient(
                radius + 0.5 * step, temperature + 0.5 * step * k1
            )
            k3 = self._material_temperature_gradient(
                radius + 0.5 * step, temperature + 0.5 * step * k2
            )
            k4 = self._material_temperature_gradient(radius + step, temperature + step * k3)
            temperature += step * (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0
            radius += step
            path.append((radius, temperature))
        return temperature, path

    def reduced_trajectory(
        self, start: float, stop: float, substeps_per_element: int = 16
    ) -> dict[str, Any]:
        if not self.radius[0] <= start < stop <= self.radius[-1]:
            raise ValueError("trajectory limits must be ordered within the mesh")
        knots = [start]
        knots.extend(float(r) for r in self.radius if start < r < stop)
        knots.append(stop)

        def integrate(substeps: int) -> tuple[np.ndarray, np.ndarray]:
            temperature = self._interpolate(self.tref, start)
            radii = [start]
            temperatures = [temperature]
            for left, right in zip(knots[:-1], knots[1:]):
                temperature, path = self._integrate_segment(
                    left, right, temperature, substeps
                )
                radii.extend(point[0] for point in path)
                temperatures.extend(point[1] for point in path)
            return np.asarray(radii), np.asarray(temperatures)

        coarse_radius, coarse_temperature = integrate(substeps_per_element)
        radius, temperature = integrate(2 * substeps_per_element)
        coarse_final = float(coarse_temperature[-1])
        integration_error_k = abs(float(temperature[-1]) - coarse_final) * self.delta_temperature_k
        reference = np.interp(radius, self.radius, self.tref)
        delta_k = (temperature - reference) * self.delta_temperature_k
        rho = np.interp(radius, self.radius, self.rho)
        velocity = 1.0 / (radius * radius * rho)
        dt = np.diff(radius) / (0.5 * (velocity[:-1] + velocity[1:]))
        time_weighted_rms = float(
            np.sqrt(np.sum(0.5 * (delta_k[:-1] ** 2 + delta_k[1:] ** 2) * dt) / np.sum(dt))
        )
        finite_difference_gradient = np.diff(reference) / np.diff(radius)
        equation_gradient = np.asarray(
            [
                self._material_temperature_gradient(float(r), float(t))
                for r, t in zip(
                    0.5 * (radius[:-1] + radius[1:]),
                    0.5 * (reference[:-1] + reference[1:]),
                )
            ]
        )
        reference_residual = finite_difference_gradient - equation_gradient
        return {
            "implementation_status": IMPLEMENTED,
            "execution_status": RUN,
            "evidence_scope": "REDUCED LAGRANGIAN ODE; no FE/SUPG assembly",
            "start_radius": start,
            "stop_radius": stop,
            "elapsed_nondimensional_time": float(np.sum(dt)),
            "samples": int(len(radius)),
            "initial_temperature_K": float(
                self._interpolate(self.tref, start) * self.delta_temperature_k
                + self.surface_temperature * self.delta_temperature_k
            ),
            "final_temperature_K": float(
                temperature[-1] * self.delta_temperature_k
                + self.surface_temperature * self.delta_temperature_k
            ),
            "maximum_deltaT_K": float(np.max(np.abs(delta_k))),
            "final_deltaT_K": float(delta_k[-1]),
            "time_weighted_RMS_deltaT_K": time_weighted_rms,
            "integration_error_estimate_K": integration_error_k,
            "reference_equation_residual_RMS": float(
                np.sqrt(np.mean(reference_residual * reference_residual))
            ),
            "state_changed_K": float(
                abs(temperature[-1] - temperature[0]) * self.delta_temperature_k
            ),
        }

    def reduced_energy_residual(
        self,
        temperature: np.ndarray | None = None,
        temperature_rate: np.ndarray | None = None,
        *,
        phases: Iterable[int] = (),
        velocity_amplitude: float = 1.0,
    ) -> tuple[EnergyTerms, np.ndarray]:
        """Cell-centred D_phase-multiplied material residual.

        This is an element-interior algebra check.  It omits FE quadrature,
        discontinuous-D interfaces, lumped mass, boundary stripping, and SUPG.
        """
        total = self.tref if temperature is None else np.asarray(temperature)
        rate = (
            np.zeros(len(self.radius) - 1)
            if temperature_rate is None
            else np.asarray(temperature_rate)
        )
        if rate.shape != (len(self.radius) - 1,):
            raise ValueError("temperature_rate must be cell centred")
        average = lambda value: 0.5 * (value[:-1] + value[1:])
        rho = average(self.rho)
        gravity = average(self.gravity)
        alpha = average(self.alpha)
        cp = average(self.cp)
        theta = average(total) + self.surface_temperature
        velocity = velocity_amplitude / (average(self.radius) ** 2 * rho)
        gradient = np.diff(total) / np.diff(self.radius)
        phase_work = np.zeros_like(rho)
        inertia = np.ones_like(rho)
        theta_nodes = total + self.surface_temperature
        nodal_velocity = velocity_amplitude / (self.radius * self.radius * self.rho)
        phase_indices = tuple(phases)
        if phase_indices:
            self._load_phase_parameters()
        for index in phase_indices:
            fraction = self.phase_fraction(index, total)
            phase_temperature = average(
                fraction * (1.0 - fraction) * theta_nodes
            )
            phase_temperature_velocity = average(
                fraction * (1.0 - fraction) * theta_nodes * nodal_velocity
            )
            coefficient = (
                2.0
                * self.dissipation_number
                * self.phase_ra[index]
                / (self.rayleigh * self.phase_width[index])
            )
            phase_work += (
                rho
                * alpha
                * gravity
                * phase_temperature_velocity
                * self.phase_clapeyron[index]
                * coefficient
            )
            inertia += (
                phase_temperature
                * self.phase_clapeyron[index] ** 2
                * coefficient
            )
        storage = inertia * rho * cp * rate
        advection = inertia * rho * cp * velocity * gradient
        adiabatic = (
            self.dissipation_number * rho * alpha * gravity * theta * velocity
        )
        zeros = np.zeros_like(rho)
        terms = EnergyTerms(
            R_storage=storage,
            R_advection=advection,
            R_adiabatic=adiabatic,
            R_phase=phase_work,
            R_diffusion=zeros.copy(),
            R_viscous=zeros.copy(),
            R_internal=zeros.copy(),
            R_boundary=zeros.copy(),
            R_stabilization=zeros.copy(),
        )
        return terms, inertia

    def reference_trajectory_diagnostics(
        self, *, include_phases: bool = True
    ) -> dict[str, Any]:
        terms, inertia = self.reduced_energy_residual()
        diagnostics = self.temperature_diagnostics(include_phases=include_phases)
        diagnostics["energy_contributions"] = {
            name: np.asarray(value).tolist() for name, value in asdict(terms).items()
        }
        diagnostics["energy_contributions"]["R_total"] = terms.total.tolist()
        diagnostics["D_phase"] = inertia.tolist()
        diagnostics["scope"] = "REDUCED CELL-CENTRED, NO FE/SUPG"
        return diagnostics

    def smooth_closure(self) -> list[dict[str, Any]]:
        """Probe the four phase-free branches of the serialized Tref profile.

        Phase terms are disabled in the ODE.  The trajectory endpoints are
        selected from the actual reference-state nodes between nominal phase
        depths; no phase law, width, Clapeyron slope, or phase amplitude is
        sampled before the smooth gate has passed.
        """
        max_depth_km = float(np.max(self.depth) * self.radius_m / 1000.0)
        depth_bounds = (0.0, *BOUNDARIES_KM, max_depth_km)
        nodal_depth_km = self.depth * self.radius_m / 1000.0
        output = []
        for index, (shallow, deep) in enumerate(zip(depth_bounds[:-1], depth_bounds[1:])):
            branch_indices = np.flatnonzero(
                (nodal_depth_km > shallow) & (nodal_depth_km < deep)
            )
            if branch_indices.size < 2:
                raise ValueError(f"smooth region {index + 1} has fewer than two nodes")
            if not np.all(np.diff(branch_indices) == 1):
                raise ValueError(f"smooth region {index + 1} nodes are not contiguous")
            element_start = int(branch_indices[0])
            element_stop = int(branch_indices[-1])
            selected_boundaries = BOUNDARIES_KM[
                (BOUNDARIES_KM >= nodal_depth_km[element_stop])
                & (BOUNDARIES_KM <= nodal_depth_km[element_start])
            ]
            if selected_boundaries.size:
                raise ValueError(
                    f"smooth region {index + 1} brackets phase boundary "
                    f"{selected_boundaries.tolist()}"
                )
            start = float(self.radius[element_start])
            stop = float(self.radius[element_stop])
            start_depth = (1.0 - start) * self.radius_m / 1000.0
            stop_depth = (1.0 - stop) * self.radius_m / 1000.0
            result = self.reduced_trajectory(start, stop)
            terms, _ = self.reduced_energy_residual()
            element_slice = slice(element_start, element_stop)
            residual = terms.total[element_slice]
            widths = np.diff(self.radius)[element_slice]
            component_norms = {
                name: _norms(np.asarray(value)[element_slice])
                for name, value in asdict(terms).items()
            }
            largest_component = max(
                component_norms,
                key=lambda name: component_norms[name]["RMS"],
            )
            result.update(
                {
                    "region": index + 1,
                    "depth_range_km": [shallow, deep],
                    "trajectory_depth_range_km": [start_depth, stop_depth],
                    "mesh_node_indices_one_based": [
                        element_start + 1,
                        element_stop + 1,
                    ],
                    "mesh_element_indices_one_based": [
                        element_start + 1,
                        element_stop,
                    ],
                    "selected_phase_boundaries_km": selected_boundaries.tolist(),
                    "trajectory_selection": (
                        "COMPLETE CONTIGUOUS SAME-BRANCH MESH ELEMENTS"
                    ),
                    "integration_tolerance_K": TRAJECTORY_INTEGRATION_TOLERANCE_K,
                    "closure_tolerance_K": TRAJECTORY_CLOSURE_TOLERANCE_K,
                    "tolerance_provenance": (
                        "ORACLE-DEFINED; request supplied no numeric trajectory "
                        "closure or integration threshold"
                    ),
                    "integrated_reduced_residual": float(
                        np.sum(residual * widths)
                    ),
                    "integrated_reduced_residual_definition": (
                        "sum_e R_e*Delta_r_e; nondimensional radial line sum, "
                        "not FE-volume, physical-power, or time integrated"
                    ),
                    "reduced_residual_norms": _norms(residual),
                    "reduced_residual_relative_RMS": _relative_rms(
                        residual, terms.R_advection[element_slice]
                    ),
                    "largest_reduced_component_by_unweighted_RMS": largest_component,
                    "largest_component_interpretation": (
                        "magnitude ranking only; cancellation prevents root-cause inference"
                    ),
                    "reduced_component_norms": component_norms,
                }
            )
            integration_ok = (
                result["integration_error_estimate_K"]
                <= TRAJECTORY_INTEGRATION_TOLERANCE_K
            )
            result["result"] = (
                PASS
                if integration_ok
                and result["maximum_deltaT_K"] <= TRAJECTORY_CLOSURE_TOLERANCE_K
                and result["time_weighted_RMS_deltaT_K"]
                <= TRAJECTORY_CLOSURE_TOLERANCE_K
                else FAIL
            )
            output.append(result)
        return output

    def buoyancy_projection_audit(self) -> dict[str, Any]:
        horizontal = np.linspace(0.0, 2.0 * math.pi, 48, endpoint=False)
        perturbation = 0.03 * np.sin(horizontal)[:, None] * np.cos(
            np.linspace(0.0, 3.0 * math.pi, len(self.radius))
        )[None, :]
        total = self.tref[None, :] + perturbation
        coefficient = self.rho[None, :] * self.alpha[None, :]
        results = []
        for partitions in (1, 2, 3, 4, 6, 8):
            absolute = horizontal_projection(coefficient * total, partitions)
            anomaly = horizontal_projection(
                coefficient * (total - self.tref[None, :]), partitions
            )
            results.append(
                {
                    "partitions": partitions,
                    "difference": _norms(absolute - anomaly),
                    "evidence_scope": (
                        "REDUCED SERIAL ARITHMETIC"
                        if partitions == 1
                        else "SIMULATED REDUCTION ORDER"
                    ),
                }
            )
        return {
            "implementation_status": IMPLEMENTED,
            "execution_status": RUN,
            "projection": "UNWEIGHTED; not production FE area projection",
            "results": results,
            "production_MPI": {
                "implementation_status": NOT_IMPLEMENTED,
                "execution_status": NOT_RUN,
                "result": None,
                "measurements": None,
                "required_world_ranks": self.mpi_ranks,
                "required_horizontal_comm_size": self.horizontal_comm_size,
            },
        }

    def zero_reference_buoyancy(self) -> dict[str, Any]:
        columns = 48
        coefficient_thermal = horizontal_projection(
            np.broadcast_to(
                self.rho * self.alpha * self.tref,
                (columns, len(self.radius)),
            )
        )
        scaled_reference = np.broadcast_to(
            self.rayleigh * self.rho * self.alpha * self.tref * self.gravity,
            (columns, len(self.radius)),
        )
        scaled_thermal = horizontal_projection(
            np.broadcast_to(
                self.rayleigh * self.rho * self.alpha * self.tref * self.gravity,
                (columns, len(self.radius)),
            )
        )
        phase_fields = []
        for index in range(3):
            fraction = self.phase_fraction(index, self.tref)
            reference = self.phase_fraction(index, self.tref)
            phase_fields.append(
                -self.phase_ra[index]
                * (fraction - reference)[None, :]
                * self.gravity[None, :]
            )
        phase = sum(phase_fields, np.zeros_like(scaled_thermal))
        nonzero_ratios = self.buoyancy_ratios[self.buoyancy_ratios != 0.0]
        if nonzero_ratios.size == 0:
            raise ValueError("zero-reference composition probe needs a nonzero ratio")
        composition_reference = 0.2 + 0.1 * np.cos(np.pi * self.depth / np.max(self.depth))
        composition_unscaled = np.broadcast_to(
            -nonzero_ratios[0] * composition_reference,
            (columns, len(self.radius)),
        )
        coefficient_composition = horizontal_projection(composition_unscaled)
        composition = horizontal_projection(
            self.rayleigh * composition_unscaled * self.gravity[None, :]
        )
        combined = scaled_thermal + phase + composition
        coefficient_norms = _norms(coefficient_thermal)
        coefficient_phase_norms = _norms(phase)
        coefficient_composition_norms = _norms(coefficient_composition)
        scaled_norms = _norms(scaled_thermal)
        phase_norms = _norms(phase)
        composition_norms = _norms(composition)
        combined_norms = _norms(combined)
        scale_norms = _norms(scaled_reference)
        normalized = _relative_norms(scaled_norms, scale_norms)
        coefficient_fields = (
            coefficient_norms,
            coefficient_phase_norms,
            coefficient_composition_norms,
        )
        scaled_fields = (scaled_norms, phase_norms, composition_norms, combined_norms)
        return {
            "implementation_status": IMPLEMENTED,
            "execution_status": RUN,
            "evidence_scope": "REDUCED SERIAL ARITHMETIC; no FE area projection",
            "serial": {
                "coefficient_space": {
                    "thermal": coefficient_norms,
                    "phase": coefficient_phase_norms,
                    "composition": coefficient_composition_norms,
                    "composition_unprojected_signal_Linf": float(
                        np.max(np.abs(composition_unscaled))
                    ),
                    "tolerance": SERIAL_BUOYANCY_TOLERANCE,
                    "result": (
                        PASS
                        if all(
                            norms["L2"] <= SERIAL_BUOYANCY_TOLERANCE
                            and norms["Linf"] <= SERIAL_BUOYANCY_TOLERANCE
                            for norms in coefficient_fields
                        )
                        else FAIL
                    ),
                },
                "production_scaled": {
                    "thermal_raw": scaled_norms,
                    "phase_raw": phase_norms,
                    "composition_raw": composition_norms,
                    "combined_raw": combined_norms,
                    "reference_force_scale": scale_norms,
                    "thermal_normalized": normalized,
                    "raw_absolute_tolerance": SERIAL_BUOYANCY_TOLERANCE,
                    "raw_absolute_result": (
                        PASS
                        if all(
                            norms["L2"] <= SERIAL_BUOYANCY_TOLERANCE
                            and norms["Linf"] <= SERIAL_BUOYANCY_TOLERANCE
                            for norms in scaled_fields
                        )
                        else FAIL
                    ),
                    "normalized_tolerance": SERIAL_BUOYANCY_TOLERANCE,
                    "normalized_result": (
                        PASS
                        if normalized["L2"] <= SERIAL_BUOYANCY_TOLERANCE
                        and normalized["Linf"] <= SERIAL_BUOYANCY_TOLERANCE
                        else FAIL
                    ),
                },
            },
            "production_MPI": {
                "implementation_status": NOT_IMPLEMENTED,
                "execution_status": NOT_RUN,
                "result": None,
                "measurements": None,
                "tolerance": MPI_BUOYANCY_TOLERANCE,
                "normalization": "required for fully scaled production field",
            },
        }

    def beta_audit(self) -> dict[str, Any]:
        interval = self.intervals[:, 3]
        nodal_average = 0.5 * (self.refstate[:-1, 5] + self.refstate[1:, 5])
        depth = 1.0 - self.radius
        boundaries = np.asarray((410.0, 520.0, 660.0)) * 1000.0 / self.radius_m
        segment = np.searchsorted(boundaries, depth, side="right")
        branch_beta_lines = [
            np.polyfit(depth[segment == sid], self.refstate[segment == sid, 5], 1)
            for sid in range(4)
        ]
        expected_rows = []
        for index in range(len(interval)):
            shallow = float(depth[index + 1])
            deep = float(depth[index])
            cuts = [shallow]
            cuts.extend(float(value) for value in boundaries if shallow < value < deep)
            cuts.append(deep)
            integral = 0.0
            for start, stop in zip(cuts[:-1], cuts[1:]):
                sid = int(np.searchsorted(boundaries, 0.5 * (start + stop), side="right"))
                line = branch_beta_lines[sid]
                integral += 0.5 * (
                    np.polyval(line, start) + np.polyval(line, stop)
                ) * (stop - start)
            expected_rows.append(integral / (deep - shallow))
        expected = np.asarray(expected_rows)
        density_secant = -np.diff(np.log(self.rho)) / np.diff(self.radius)
        material = (GLOBAL_ROOT / "lib" / "Material_properties.c").read_text(
            encoding="utf-8"
        )
        element = (GLOBAL_ROOT / "lib" / "Element_calculations.c").read_text(
            encoding="utf-8"
        )
        vsolver = self.cfg["CitcomS.solver.vsolver"]
        param = self.cfg["CitcomS.solver.param"]
        selected = vsolver.get("ala_beta_element_source")
        interval_name = param.get("ala_beta_interval_file")
        inventory = [
            {
                "field": "beta_ala",
                "role": "nodal refstate column 6",
                "selected": False,
                "runtime_consumer": "candidate construction/diagnostics",
            },
            {
                "field": "ala_beta_supplied",
                "role": "endpoint-average candidate",
                "selected": False,
                "runtime_consumer": "causal diagnostic only",
            },
            {
                "field": "ala_beta_density",
                "role": "log-density-secant candidate",
                "selected": False,
                "runtime_consumer": "causal diagnostic only",
            },
            {
                "field": "ala_beta_interval",
                "role": "serialized interval candidate",
                "selected": selected == "interval",
                "runtime_consumer": "selected into ala_beta",
            },
            {
                "field": "ala_beta",
                "role": "single selected runtime authority",
                "selected": True,
                "runtime_consumer": "continuity, transpose, AL, diagnostics",
            },
            {
                "field": "gamma_eff",
                "role": "thermodynamic closure diagnostic, not beta",
                "selected": False,
                "runtime_consumer": "reference validation only",
            },
        ]
        contracts = {
            "cfg_selects_interval": selected == "interval",
            "cfg_selects_file_reference": param.getint("reference_state") == 0,
            "cfg_interval_file_matches_input": interval_name == self.inputs.interval.name,
            "interval_reader": "read_ala_beta_intervals(E);" in material,
            "selected_runtime_assignment": "E->refstate.ala_beta[nz] = beta;" in material,
            "continuity_uses_selected_beta": "E->refstate.ala_beta[fine_nz] * dr" in element,
            "transpose_uses_same_element_matrix": (
                "ALA_COMBINED_PRESSURE_COEFFICIENT(" in element
                and "E->elt_c[lev][m][e].c[p][0]" in element
            ),
        }
        compiled_sources = {
            str(path.relative_to(GLOBAL_ROOT)): path.read_text(encoding="utf-8")
            for path in _compiled_library_sources()
        }
        source_scan = scan_beta_source_occurrences(compiled_sources)
        declaration_inventory = ref_state_beta_gamma_declarations(
            compiled_sources["lib/global_defs.h"]
        )
        declarations_match = (
            declaration_inventory == BETA_REF_STATE_DECLARATION_CONTRACT
        )
        contracts["frozen_occurrence_multiset_matches"] = source_scan[
            "contract_matches"
        ]
        contracts["refstate_beta_gamma_declarations_frozen"] = declarations_match
        return {
            "implementation_status": IMPLEMENTED,
            "execution_status": RUN,
            "selected_cfg_source": selected,
            "serialized_vs_regenerated": _norms(interval - expected),
            "supplied_average_vs_regenerated": _norms(nodal_average - expected),
            "density_secant_vs_regenerated": _norms(density_secant - expected),
            "inventory": inventory,
            "source_contracts": contracts,
            "source_scan": {
                "inventory_source": "lib/Makefile.am compiled sources",
                **source_scan,
                "declaration_inventory": list(declaration_inventory),
                "expected_declaration_inventory": list(
                    BETA_REF_STATE_DECLARATION_CONTRACT
                ),
                "declaration_inventory_matches": declarations_match,
                "known_field_scope": (
                    "All beta/Gamma-like REF_STATE member declarations are frozen; "
                    "the scan does not prove absence of locally regenerated beta "
                    "that never uses a REF_STATE beta/Gamma member"
                ),
            },
            "unexpected_runtime_authorities": source_scan[
                "unexpected_occurrences"
            ],
            "missing_expected_occurrences": source_scan["missing_occurrences"],
        }

    @staticmethod
    def blocked_experiment(reason: str) -> dict[str, Any]:
        return {
            "implementation_status": NOT_IMPLEMENTED,
            "execution_status": BLOCKED,
            "result": None,
            "measurements": None,
            "reason": reason,
        }

    def source_contracts(self) -> dict[str, bool]:
        energy = (GLOBAL_ROOT / "lib" / "Advection_diffusion.c").read_text(
            encoding="utf-8"
        )
        phase = (GLOBAL_ROOT / "lib" / "Phase_change.c").read_text(encoding="utf-8")
        buoyancy = (GLOBAL_ROOT / "lib" / "Pan_problem_misc_functions.c").read_text(
            encoding="utf-8"
        )
        return {
            "total_temperature_is_energy_state": "pg_solver(E,E->T,E->Tdot" in energy,
            "phase_uses_total_temperature": (
                "phase_change_state(phase" in phase and "E->T[m][i]" in phase
            ),
            "phase_reference_uses_Tref": (
                "phase_change_state(phase" in phase and "E->refstate.Tref[nz]" in phase
            ),
            "phase_buoyancy_is_anomaly": "phase->Ra * (B[m][i] - Xref)" in phase,
            "phase_energy_uses_single_entropy": "phase->entropy_jump" in energy,
            "phase_energy_uses_current_temperature": "tgp[i]" in energy,
            "phase_work_uses_complete_proxy_derivative": (
                "(depth - phase->depth) * d_rho_g_dr" in phase
            ),
            "obsolete_phase_operator_is_inactive": (
                "static void latent_heating" not in energy
                and "E->phase_B[phase_index]" not in energy
            ),
            "SUPG_active_with_zero_diffusion": "unorm>0.000001" in energy,
            "horizontal_projection_is_applied": "remove_horiz_ave2(E,buoy);" in buoyancy,
            "rheo_dat_can_mutate_temperature": (
                'fp=fopen("rheo.dat","r");' in energy
                and "E->T[m][i] +=" in energy
            ),
            "no_explicit_DXDt_source": not bool(re.search(r"(?:DB|DX)\s*/\s*Dt", energy)),
        }

    def run(self) -> dict[str, Any]:
        # Required prerequisite audit: source text only, with no phase inputs loaded.
        source_contracts = self.source_contracts()
        beta = self.beta_audit()
        projection = self.buoyancy_projection_audit()
        smooth = self.smooth_closure()
        reduced_smooth_pass = all(item["result"] == PASS for item in smooth)
        full_fe_smooth_pass = None
        phase_reason = (
            "Phase crossing requires smooth reduced and production-FE closure; "
            f"reduced={PASS if reduced_smooth_pass else FAIL}, production-FE={NOT_IMPLEMENTED}"
        )
        phase_closure = {
            name: self.blocked_experiment(phase_reason) for name in PHASE_NAMES
        }
        count_once = {
            name: self.blocked_experiment(
                "Requires controlled evolved phase/temperature toggles after smooth closure"
            )
            for name in ("experiment_A", "experiment_B", "experiment_C")
        }
        forward_reverse = {
            name: self.blocked_experiment(
                "Requires a forward/reverse evolved trajectory after smooth closure"
            )
            for name in PHASE_NAMES
        }
        phase_ledger = {
            "implementation_status": NOT_IMPLEMENTED,
            "execution_status": BLOCKED,
            "result": None,
            "values": None,
            "reason": phase_reason,
        }
        unavailable_fe = {
            "implementation_status": NOT_IMPLEMENTED,
            "execution_status": NOT_RUN,
            "result": None,
            "measurements": None,
        }
        if reduced_smooth_pass and full_fe_smooth_pass:
            ownership = self.reference_trajectory_diagnostics(include_phases=True)
            zero = self.zero_reference_buoyancy()
            coefficient_zero_gate = zero["serial"]["coefficient_space"]["result"]
            raw_zero_gate = zero["serial"]["production_scaled"]["raw_absolute_result"]
            normalized_zero_gate = zero["serial"]["production_scaled"]["normalized_result"]
        else:
            ownership = self.reference_trajectory_diagnostics(include_phases=False)
            zero = {
                "implementation_status": NOT_IMPLEMENTED,
                "execution_status": BLOCKED,
                "result": None,
                "measurements": None,
                "reason": phase_reason,
            }
            coefficient_zero_gate = BLOCKED
            raw_zero_gate = BLOCKED
            normalized_zero_gate = BLOCKED
        return {
            "oracle_scope": "TEST-ONLY reduced/static evidence framework",
            "production_files_modified": False,
            "manifest": self.input_manifest(),
            "temperature_ownership": ownership,
            "manufactured_mass_flux": self.mass_flux_constancy(),
            "smooth_closure": smooth,
            "phase_closure": phase_closure,
            "phase_energy_ledger": phase_ledger,
            "count_once": count_once,
            "thermal_buoyancy_equivalence": projection,
            "zero_reference_buoyancy": zero,
            "beta_consistency": beta,
            "forward_reverse": forward_reverse,
            "source_contracts": source_contracts,
            "source_contracts_scope": (
                "PREREQUISITE STATIC ENERGY MAP; no phase parameter loading, "
                "phase-law evaluation, or phase experiment"
            ),
            "unavailable_production_evidence": {
                "smooth_FE_trajectory": dict(unavailable_fe),
                "full_FE_residual_terms": dict(unavailable_fe),
                "production_MPI_projection": dict(unavailable_fe),
            },
            "gates": {
                "serial_zero_reference_coefficient_space": (
                    coefficient_zero_gate
                ),
                "serial_zero_reference_raw_production_scale": (
                    raw_zero_gate
                ),
                "serial_zero_reference_normalized_production_scale": (
                    normalized_zero_gate
                ),
                "smooth_reduced_trajectory": PASS if reduced_smooth_pass else FAIL,
                "smooth_production_FE_trajectory": NOT_IMPLEMENTED,
                "phase_crossings": BLOCKED,
                "production_MPI_projection": NOT_IMPLEMENTED,
                "full_FE_residual_terms": NOT_IMPLEMENTED,
            },
            "root_cause_classification": (
                "Reduced smooth-reference trajectory does not close"
                if not reduced_smooth_pass
                else "Production FE/SUPG evidence is unavailable"
            ),
            "verdict": "UNRESOLVED",
            "recommended_next_action": (
                "Complete the test-only initialized CitcomS FE/SUPG trajectory fixture; "
                "do not modify production physics or phase/reference parameters"
            ),
            "internal_gate_state": {
                "reduced_smooth_pass": reduced_smooth_pass,
                "full_fe_smooth_pass": full_fe_smooth_pass,
            },
        }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--refstate", type=Path, default=DEFAULT_REFSTATE)
    parser.add_argument("--interval", type=Path, default=DEFAULT_INTERVAL)
    parser.add_argument("--mesh", type=Path, default=DEFAULT_MESH)
    parser.add_argument("--katsura", type=Path, default=DEFAULT_KATSURA)
    parser.add_argument("--json", action="store_true", help="emit the complete JSON report")
    parser.add_argument(
        "--allow-unresolved",
        action="store_true",
        help="emit an unresolved diagnostic report with a zero process status",
    )
    args = parser.parse_args()
    try:
        oracle = StrictALAClosureOracle(
            OracleInputs(
                args.config,
                args.refstate,
                args.interval,
                args.mesh,
                args.katsura,
            )
        )
        report = oracle.run()
    except Exception as error:
        report = {
            "verdict": "UNRESOLVED",
            "root_cause_classification": "ORACLE-EXECUTION-FAILURE",
            "error": {"type": type(error).__name__, "message": str(error)},
            "gates": {"oracle_execution": FAIL},
            "production_files_modified": False,
        }
    if args.json:
        print(json.dumps(report, indent=2, sort_keys=True))
    else:
        print(f"Stage 0 thermodynamic verdict: {report['verdict']}")
        print(f"Root cause: {report['root_cause_classification']}")
        for item in report.get("smooth_closure", []):
            print(
                f"smooth region {item['region']}: {item['result']} "
                f"max_deltaT={item['maximum_deltaT_K']:.9g} K "
                f"RMS_deltaT={item['time_weighted_RMS_deltaT_K']:.9g} K "
                f"integration_error={item['integration_error_estimate_K']:.3e} K"
            )
        if "error" in report:
            print(
                f"oracle error: {report['error']['type']}: "
                f"{report['error']['message']}"
            )
        else:
            print("phase crossings: BLOCKED by smooth production-FE gate")
    return 0 if report["verdict"] == PASS or args.allow_unresolved else 2


if __name__ == "__main__":
    raise SystemExit(main())
