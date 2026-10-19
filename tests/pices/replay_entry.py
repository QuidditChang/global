"""Reuse the checkpoint producer's normalization in the isolated C replay.

The Pyre and standalone entries have algebraically equivalent but not bitwise
equivalent Rayleigh expressions. Copy the producer expression, including its
float radius in metres, before any phase/mesh initialization. Never patch a
checkpoint hash or the production source. Fail closed if either anchor changes.
"""
import re

VERSION = 1


def normalization(source):
    pattern = (r'(?m)^ +if\(E->data\.ref_temperature > 0\.0\)\n'
               r' +E->control\.Atemp =.*?'
               r'(?=\n\n +if\(E->data\.Ttop < 0\.0\))')
    matches = list(re.finditer(pattern, source, re.S))
    if len(matches) != 1:
        raise ValueError('Expected exactly one reference normalization block')
    return matches[0]


def adapt(instructions, properties):
    setter = properties.split('PyObject * pyCitcom_Const_set_properties(', 1)[1]
    setter = setter.split('char pyCitcom_IC_set_properties__doc__', 1)[0]
    if 'float radius;' not in setter:
        raise ValueError('Producer radius type changed')
    original = normalization(instructions)
    producer = normalization(setter).group()
    anchor = 'input_float("radius",&tmp,"6371e3.0",m);'
    if instructions.count(anchor) != 1:
        raise ValueError('Standalone radius input changed')
    between = instructions.split(anchor, 1)[1][:original.start() - instructions.index(anchor) - len(anchor)]
    if re.search(r'\btmp\s*=|&tmp', between):
        raise ValueError('Standalone radius temporary was overwritten')
    if 'radius * radius * radius' not in producer or '1e9' not in original.group():
        raise ValueError('Normalization expressions changed; review entry adapter')
    block = re.sub(r'\bradius\b', 'tmp', producer)
    return instructions[:original.start()] + block + instructions[original.end():]
