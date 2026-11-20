#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
"""Add a public LossyMetal control to the standard benchmark runner.

python3 benchmark/lossymetal_cpu.py run --cases ideal,lossymetal --size XS ...
python3 benchmark/lossymetal_cpu.py compare <result.json>

The LossyMetal case is XS only. Its sinusoidal source keeps it running until
bench.py requests ABORT. Uses the same Python bindings as benchmark/bench.py.
"""
from pathlib import Path
import xml.etree.ElementTree as ET
import bench


class LossyMetalControl:
    def Write2XML(self, filename):
        path = Path(__file__).resolve().parents[1] / 'tests/fixtures/lossymetal-setup/time-domain.xml'
        root = ET.parse(path).getroot()
        fdtd = root.find('FDTD')
        fdtd.set('NumberOfTimesteps', '100000000')
        fdtd.set('endCriteria', '1e-30')
        pulse = fdtd.find('Excitation')
        pulse.attrib.clear()
        pulse.attrib.update(Type='1', f0='20000000000')
        ET.ElementTree(root).write(filename, encoding='utf-8', xml_declaration=True)


def lossymetal(size):
    if size != 'XS':
        raise ValueError('LossyMetal control supports XS only')
    return LossyMetalControl()


if __name__ == '__main__':
    bench.cases.CASES['lossymetal'] = lossymetal
    bench.cases.CASES_VERSION = 1003
    bench.main()
