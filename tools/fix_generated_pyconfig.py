#!/usr/bin/env python3
"""Repair pyconfig emitted by an existing, pre-fix configure script.

The m4 source has the same fix. This covers HPC checkouts whose untracked
configure script predates a git pull, without requiring an Autotools upgrade.
"""
from pathlib import Path
import sys

OLD = '''keys = makefile_vars.keys()
for key in keys:
    makefile_vars[key] = expand_makefile_vars(makefile_vars[key], makefile_vars)
'''
NEW = '# parse_makefile already expands variables and preserves integer values.\n'


def repair(path):
    source = path.read_text()
    if OLD in source:
        path.write_text(source.replace(OLD, NEW))
    elif NEW not in source:
        raise ValueError("unrecognized pyconfig template: " + str(path))


if __name__ == '__main__':
    repair(Path(sys.argv[1]))
