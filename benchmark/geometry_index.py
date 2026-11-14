#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Krzysztof Rosinski
"""Add a public polygon-heavy XS setup case to the unchanged benchmark runner.

Example: geometry_index.py run --openems /path/to/master --label master
    --openems /path/to/indexed-wrapper --label indexed --cases polygons
    --size XS --engine basic --threads 2 --duration .5 --repeat 3

The indexed wrapper sets OPENEMS_GEOMETRY_INDEX=16,16,1. No index is enabled
implicitly by this script. The extra case is a setup workload, not an antenna.
"""
from pathlib import Path
import runpy

import numpy as np

import cases


def polygons(size):
    fdtd = cases.inhomogeneous(size)
    csx = fdtd.GetCSX()
    extent = float(csx.GetGrid().GetLines('x')[-1])
    metal = csx.AddMetal('public_polygon_workload')
    pitch = 1.4 * extent / 32
    angle = np.linspace(0, 2 * np.pi, 33)
    for x in range(32):
        for y in range(32):
            centre = np.array([-.7 * extent + pitch * (x + .5), -.7 * extent + pitch * (y + .5)])
            points = centre[:, None] + .35 * pitch * np.array([np.cos(angle), np.sin(angle)])
            metal.AddLinPoly(points=points, norm_dir=2, elevation=.2 * extent,
                             length=.25, priority=20 + (x + y) % 7)
    return fdtd


if __name__ == '__main__':
    cases.CASES['polygons'] = polygons
    cases.CASES_VERSION = '5+public-polygons-1024'
    runpy.run_path(str(Path(__file__).with_name('bench.py')), run_name='__main__')
