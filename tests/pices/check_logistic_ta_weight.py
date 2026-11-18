"""Check the compiled TA weight against independent Decimal quadrature.

Run with ``python3 tests/pices/check_logistic_ta_weight.py``.  Only a C compiler
and the Python standard library are needed; all build files are temporary.
The reference integrates the prescribed logistic slope, rather than repeating
the production softplus antiderivative or its cancellation-avoidance series.
"""

import ctypes
from decimal import Decimal, localcontext
from functools import lru_cache
import math
import os
from pathlib import Path
import shlex
import subprocess
import tempfile


ROOT = Path(__file__).resolve().parents[2]
D = Decimal
H_KM = 127.42
CENTER = 80.0 / H_KM
WIDTH = 4.0 / H_KM


@lru_cache(maxsize=None)
def gauss_rule(order):
    """Generate Gauss-Legendre nodes/weights at the reference precision."""
    rule = []
    for index in range(1, order // 2 + 1):
        root = D.from_float(math.cos(math.pi * (index - 0.25) / (order + 0.5)))
        for _ in range(30):
            previous, value = D(1), root
            for degree in range(2, order + 1):
                previous, value = value, (
                    (2 * degree - 1) * root * value - (degree - 1) * previous
                ) / degree
            derivative = order * (root * value - previous) / (root * root - 1)
            correction = value / derivative
            root -= correction
            if abs(correction) < D("1e-65"):
                break
        else:
            raise AssertionError("reference quadrature root did not converge")
        weight = 2 / ((1 - root * root) * derivative * derivative)
        rule.extend(((-root, weight), (root, weight)))
    return tuple(rule)


def logistic_decimal(position, center, width):
    argument = (center - position) / width
    # Exponentiate only nonpositive arguments, including very narrow profiles.
    if argument >= 0:
        return 1 / (1 + (-argument).exp())
    exponential = argument.exp()
    return exponential / (1 + exponential)


def quadrature(function, left, right, level=0):
    # A monotone integrand below this bound cannot affect a binary64 result
    # for any tested normalization.  This avoids resolving exponentially tiny
    # tails hundreds of thousands of transition widths from their center.
    if abs(function(left)) < D("1e-400"):
        return D(0)
    midpoint, radius = (left + right) / 2, (right - left) / 2

    def estimate(order):
        return radius * sum(
            weight * function(midpoint + radius * node)
            for node, weight in gauss_rule(order)
        )

    coarse, fine = estimate(16), estimate(32)
    if fine == coarse or abs(fine - coarse) <= abs(fine) * D("2e-23"):
        return fine
    if level >= 24:
        raise AssertionError("independent reference quadrature did not converge")
    return quadrature(function, left, midpoint, level + 1) + quadrature(
        function, midpoint, right, level + 1
    )


@lru_cache(maxsize=None)
def slope_integral(fraction, center, width):
    with localcontext() as context:
        context.prec = 70
        left, c, s = map(D.from_float, (fraction, center, width))
        bottom = logistic_decimal(D(1), c, s)

        def slope(position):
            return logistic_decimal(position, c, s) - bottom

        # Resolve a narrow transition explicitly so quadrature cannot overlook
        # it between samples.  These are integration breakpoints, not a formula
        # for the reference result.
        breaks = [left, D(1)]
        breaks.extend(
            point for factor in (-40, -12, -4, 0, 4, 12, 40, 80, 160, 320, 640)
            if left < (point := c + factor * s) < 1
        )
        breaks = sorted(set(breaks))
        return sum(quadrature(slope, a, b) for a, b in zip(breaks, breaks[1:]))


def reference(fraction, center, width):
    if fraction <= 0:
        return 1.0
    if fraction >= 1:
        return 0.0
    with localcontext() as context:
        context.prec = 70
        return float(slope_integral(fraction, center, width)
                     / slope_integral(0.0, center, width))


def close(actual, expected, label, relative=2e-10, absolute=2e-14):
    assert math.isfinite(actual), (label, "nonfinite", actual)
    assert math.isclose(actual, expected, rel_tol=relative, abs_tol=absolute), (
        label, actual, expected
    )


def check_weight(weight):
    # The physical proposal, including the very weak but nonzero deep tail.
    for depth in (0, 14.564106, 25, 50, 60, 70, 75, 80, 90, 100,
                  101.948742, 110, 116.512848, 120, 125, H_KM):
        fraction = depth / H_KM
        close(weight(fraction, CENTER, WIDTH), reference(fraction, CENTER, WIDTH),
              ("proposal depth km", depth))

    # Independent published proposal checkpoints, useful when reading failures.
    close(weight(50 / H_KM, CENTER, WIDTH), 0.3750247935, "50 km", absolute=1e-10)
    close(weight(100 / H_KM, CENTER, WIDTH), 0.0003329816437,
          "100 km", absolute=1e-13)

    # Endpoints are exact; the preceding representable interior points must
    # remain positive.  Absolute tolerances would hide cancellation to zero.
    for fraction in (-1.0, 0.0):
        assert weight(fraction, CENTER, WIDTH) == 1.0
    for fraction in (1.0, 2.0):
        assert weight(fraction, CENTER, WIDTH) == 0.0
    for distance in (1e-3, 1e-5, 1e-7, 1e-9, 1e-12, 1e-15):
        fraction = 1.0 - distance
        actual = weight(fraction, CENTER, WIDTH)
        assert actual > 0.0, ("premature zero near bottom", fraction)
        close(actual, reference(fraction, CENTER, WIDTH),
              ("bottom tail", fraction), relative=2e-8, absolute=0.0)
    fraction = math.nextafter(1.0, 0.0)
    actual = weight(fraction, CENTER, WIDTH)
    assert actual > 0.0, "last representable interior point lost its tail"
    close(actual, reference(fraction, CENTER, WIDTH), "last interior point",
          relative=2e-8, absolute=0.0)

    # Profiles with a narrow transition, with a center close to the surface or
    # bottom, and with a wide transition exercise different numerical scales.
    profiles = ((CENTER, WIDTH),
                (ctypes.c_float(CENTER).value, ctypes.c_float(WIDTH).value),
                (0.6, 1e-6), (0.6, 1e-3),
                (0.02, 0.01), (0.98, 0.01), (0.6, 10.0), (0.6, 1e8))
    for center, width in profiles:
        points = {0.01, 0.25, 0.5, 0.9, 0.99, center}
        points.update(center + factor * width for factor in (-4, -1, 1, 4)
                      if 0 < center + factor * width < 1)
        for fraction in sorted(points):
            close(weight(fraction, center, width), reference(fraction, center, width),
                  ("profile", fraction, center, width), relative=3e-8)
        values = [weight(index / 4096, center, width) for index in range(4097)]
        assert all(math.isfinite(value) and 0 <= value <= 1 for value in values)
        differences = [b - a for a, b in zip(values, values[1:])]
        assert max(differences) <= 2e-13, ("weight increases", center, width)
        assert all(b >= a - 2e-12 for a, b in zip(differences, differences[1:])), (
            "slope reverses instead of flattening", center, width
        )

    # Check the defining slope equation, not an implementation-specific series.
    denominator = float(slope_integral(0.0, CENTER, WIDTH))
    bottom = 1 / (1 + math.exp((1 - CENTER) / WIDTH))
    for fraction in (0.1, 0.4, 0.6, 0.7, 0.8, 0.9, 0.99):
        half_step = 2e-6
        measured = (weight(fraction + half_step, CENTER, WIDTH)
                    - weight(fraction - half_step, CENTER, WIDTH)) / (2 * half_step)
        logistic = 1 / (1 + math.exp((fraction - CENTER) / WIDTH))
        close(measured, -(logistic - bottom) / denominator,
              ("slope equation", fraction), relative=1e-6, absolute=2e-11)

    # Center and width scale with the local H, including the protected trench
    # column; a change of length units cannot change the prescribed weight.
    for fraction in (0.1, 0.5, 0.8, 0.95):
        expected = weight(fraction, CENTER, WIDTH)
        for thickness in (H_KM, 0.003 * H_KM, H_KM * 1000):
            close(weight((fraction * thickness) / thickness,
                         (CENTER * thickness) / thickness,
                         (WIDTH * thickness) / thickness), expected,
                  ("length scaling", fraction, thickness))


def main():
    with tempfile.TemporaryDirectory(prefix="citcoms-logistic-ta-") as temporary:
        folder = Path(temporary)
        source, library = folder / "wrapper.c", folder / "weight.so"
        source.write_text(
            '#include "lith_age_weight.h"\n'
            'double evaluate(double x, double c, double s) {\n'
            '    return lith_age_logistic_rate_weight(x, c, s);\n}\n'
        )
        compiler = shlex.split(os.environ.get("CC", "cc"))
        subprocess.run(compiler + ["-std=c99", "-O2", "-Wall", "-Wextra",
                       "-shared", "-fPIC", "-I" + str(ROOT / "lib"),
                       str(source), "-lm", "-o", str(library)], check=True)
        loaded = ctypes.CDLL(str(library))
        weight = loaded.evaluate
        weight.argtypes = [ctypes.c_double] * 3
        weight.restype = ctypes.c_double
        check_weight(weight)
    print("PASS: compiled logistic TA weight matches independent high-precision "
          "quadrature; endpoints, deep tail, monotone slope, width extremes "
          "and local-depth scaling verified")


if __name__ == "__main__":
    main()
