# SPDX-License-Identifier: GPL-3.0-or-later
"""Independent continuum reference for the scoped closed-box XY packet benchmark.

Normalized E_z or eta*H_z in a 30 mm vacuum square. No openEMS import or field
fitting. Native seed time-stagger residuals remain part of any solver comparison.
"""
from pathlib import Path
import argparse
import json
import math

import numpy as np


LENGTH = .03
UNIT = .01
EPSILON = 8.85418781762e-12
MU = 1.2566370614359173e-6
SPEED = 1 / math.sqrt(EPSILON * MU)


def require(condition, message):
    if not condition:
        raise ValueError(message)



def profile(case, x, y):
    angle = math.radians(case['angle'])
    c, s = math.cos(angle), math.sin(angle)
    longitudinal = (c * x + s * y) / UNIT + .7
    transverse = (-s * x + c * y) / UNIT - .07
    return np.exp(-(longitudinal / case['width'])**2 - (transverse / .35)**2) * \
        np.cos(2 * math.pi * longitudinal / case['wavelength'])



def receivers(case):
    angle = math.radians(case['angle'])
    c, s = math.cos(angle), math.sin(angle)
    return np.array([[UNIT * (c * a - s * .07), UNIT * (s * a + c * .07)]
                     for a in [-.6, 0, .65]])



class Modes:
    def __init__(self, magnetic, maximum, quadrature, seed):
        self.magnetic = magnetic
        self.indices = np.arange(0 if magnetic else 1, maximum + 1)
        self.k = self.indices * math.pi / LENGTH
        nodes, weights = np.polynomial.legendre.leggauss(quadrature)
        self.nodes = nodes * LENGTH / 2
        self.weights = weights * LENGTH / 2
        basis = self.basis(self.nodes)
        normalization = np.full(len(self.indices), 2 / LENGTH)
        if magnetic:
            normalization[0] = 1 / LENGTH
        values = seed(self.nodes[:, None], self.nodes[None, :])
        require(values.shape == (quadrature, quadrature) and np.isfinite(values).all(), 'Invalid seed')
        self.coefficients = (basis.T @ (self.weights[:, None] * values * self.weights[None, :]) @ basis) * \
            normalization[:, None] * normalization[None, :]
        self.omega = SPEED * np.sqrt(self.k[:, None]**2 + self.k[None, :]**2)

    def basis(self, coordinates, half_width=0):
        phases = (np.asarray(coordinates)[:, None] + LENGTH / 2) * self.k[None, :]
        values = np.cos(phases) if self.magnetic else np.sin(phases)
        # Exact mean over an interval: sin(ka)/(ka); numpy sinc uses pi*x.
        return values * np.sinc(self.k[None, :] * half_width / math.pi)

    def spatial_weights(self, points, half_width):
        a, b = self.basis(points[:, 0], half_width), self.basis(points[:, 1], half_width)
        return (a[:, :, None] * b[:, None, :]).reshape(len(points), -1)

    def evaluate(self, times, points, half_width=0):
        spatial = self.spatial_weights(points, half_width)
        values = np.empty((len(times), len(points)))
        weighted = self.coefficients.ravel()[:, None] * spatial.T
        for first in range(0, len(times), 32):
            part = np.asarray(times[first:first + 32])
            values[first:first + 32] = np.cos(part[:, None] * self.omega.ravel()[None, :]) @ weighted
        require(np.isfinite(values).all(), 'Nonfinite analytical field')
        return values



def controls():
    rows = []
    points = np.array([[-.0081, .0032], [.0019, -.0028], [.0067, .0083]])
    for magnetic in [False, True]:
        a, b = (2, 3) if magnetic else (3, 5)
        trig = np.cos if magnetic else np.sin
        kx, ky = a * math.pi / LENGTH, b * math.pi / LENGTH
        seed = lambda x, y: trig(kx * (x + LENGTH / 2)) * trig(ky * (y + LENGTH / 2))
        modes = Modes(magnetic, 8, 64, seed)
        expected = np.zeros_like(modes.coefficients)
        expected[a if magnetic else a - 1, b if magnetic else b - 1] = 1
        coefficient_error = float(np.max(np.abs(modes.coefficients - expected)))
        require(coefficient_error < 1e-12, 'Modal normalization/orthogonality failed')
        frequency = SPEED * math.sqrt(kx*kx + ky*ky)
        times = np.array([0, math.pi / (2 * frequency), math.pi / frequency, 2 * math.pi / frequency])
        exact = np.cos(frequency * times[:, None]) * seed(points[:, 0], points[:, 1])[None, :]
        point_error = float(np.max(np.abs(modes.evaluate(times, points) - exact)))
        require(point_error < 1e-12, 'Known standing-mode time/position check failed')
        half = .0008
        exact_mean = exact * np.sinc(kx * half / math.pi) * np.sinc(ky * half / math.pi)
        mean_error = float(np.max(np.abs(modes.evaluate(times, points, half) - exact_mean)))
        require(mean_error < 1e-12, 'Exact volume mean failed')
        # Independent direct Gauss integral checks the endpoint/sinc formula.
        nodes, weights = np.polynomial.legendre.leggauss(32)
        direct = np.array([np.sum(weights[:, None] * weights[None, :] *
                          seed(p[0] + half * nodes[:, None], p[1] + half * nodes[None, :])) / 4
                           for p in points])
        quadrature_error = float(np.max(np.abs(direct - exact_mean[0])))
        require(quadrature_error < 1e-12, 'Independent volume integral failed')
        if magnetic:
            wall_derivative = -modes.k * np.sin(modes.k * LENGTH) * LENGTH
            boundary_error = float(np.max(np.abs(wall_derivative)))
        else:
            boundary_error = float(np.max(np.abs(modes.basis(np.array([-LENGTH / 2, LENGTH / 2])))))
        require(boundary_error < 1e-10, 'Wrong physical boundary basis')
        rows.append(dict(magnetic=magnetic, coefficientError=coefficient_error, pointError=point_error,
                         volumeError=mean_error, independentVolumeError=quadrature_error,
                         boundaryError=boundary_error))
    constant = Modes(True, 8, 64, lambda x, y: np.ones(np.broadcast_shapes(x.shape, y.shape)))
    expected = np.zeros_like(constant.coefficients)
    expected[0, 0] = 1
    require(np.max(np.abs(constant.coefficients - expected)) < 1e-12, 'Wrong zero-mode normalization')
    require(np.max(np.abs(constant.evaluate(np.array([0, 1e-9]), points, .0008) - 1)) < 1e-12,
            'Neumann DC mode evolved incorrectly')
    rows.append(dict(constantNeumannMode=True))
    return rows



def relative_difference(a, b):
    peaks = np.max(np.abs(b), axis=0)
    require(np.all(peaks > 1e-12), 'Vanishing analytical receiver reference')
    return np.max(np.abs(a - b), axis=0) / peaks



def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input', nargs='?', type=Path, help='Packet CSV with epoch,time_e,time_h')
    parser.add_argument('--output', type=Path, help='New analytical CSV; existing files are refused')
    parser.add_argument('--angle', type=float)
    parser.add_argument('--polarization', choices=['TE', 'TM'])
    parser.add_argument('--width', type=float, default=.25)
    parser.add_argument('--wavelength', type=float, default=1.5)
    parser.add_argument('--controls', action='store_true')
    args = parser.parse_args()
    if args.controls:
        require(args.input is None and args.output is None, 'Controls do not read/write packet files')
        print(json.dumps(dict(status='ANALYTICAL_BASIS_CONTROLS_PASSED', controls=controls()), indent=2))
        return
    require(args.input is not None and args.output is not None and args.angle is not None and
            args.polarization is not None, 'Input, output, angle and polarization are required')
    require(not args.output.exists(), 'Refuse to overwrite an existing output')
    require(0 <= args.angle <= 90 and .18 <= args.width <= .35 and args.wavelength == 1.5,
            'Parameters outside the scoped diagnostic range')
    data = np.genfromtxt(args.input, delimiter=',', names=True, ndmin=1)
    require(data.dtype.names is not None and {'epoch', 'time_e', 'time_h'} <= set(data.dtype.names),
            'Missing raw epoch/E/H timestamps')
    key = 'time_h' if args.polarization == 'TM' else 'time_e'
    require(1 < len(data) <= 2048 and np.isfinite(data[key]).all() and
            np.all(np.diff(data[key]) > 0) and data[key][0] >= 0 and data[key][-1] <= 81e-12,
            'Invalid or out-of-scope observation times')
    dt = data['time_e'][1] - data['time_e'][0]
    require(dt > 0 and np.array_equal(data['epoch'], np.arange(1, len(data) + 1)) and
            np.isfinite(data['time_e']).all() and np.isfinite(data['time_h']).all() and
            np.max(np.abs(data['time_e'] / dt - data['epoch'])) < 2e-8 and
            np.max(np.abs((data['time_h'] - data['time_e']) / dt - .5)) < 2e-8,
            'Expected unchanged complete-epoch E/H half-step timestamps')
    case = dict(angle=args.angle, polarization=args.polarization, width=args.width, wavelength=args.wavelength)
    points, times = receivers(case), np.r_[0., data[key]]
    traces = {}
    for maximum, quadrature in [(128, 512), (192, 384), (192, 512)]:
        modes = Modes(args.polarization == 'TM', maximum, quadrature, lambda x, y: profile(case, x, y))
        traces[(maximum, quadrature)] = np.column_stack([modes.evaluate(times, points),
                                                        modes.evaluate(times, points, .0008)])
    final = traces[(192, 512)]
    quadrature_error = float(np.max(relative_difference(traces[(192, 384)], final)))
    modal_error = float(np.max(relative_difference(traces[(128, 512)], final)))
    require(quadrature_error <= 1e-5 and modal_error <= 1e-4, 'Analytical reference did not converge')
    header = 'epoch,time_e,time_h,' + ','.join(['analytic_point_' + str(p) for p in range(3)] +
                                             ['analytic_mean_' + str(p) for p in range(3)])
    columns = np.column_stack([data['epoch'], data['time_e'], data['time_h'], final[1:]])
    require(np.isfinite(columns).all(), 'Nonfinite analytical export')
    with args.output.open('x') as stream:
        np.savetxt(stream, columns, delimiter=',', header=header, comments='', fmt='%.17g')
    print(json.dumps(dict(status='CONVERGED_CLOSED_BOX_CONTINUUM_REFERENCE_WRITTEN',
                         quadratureChange=quadrature_error, modalChange=modal_error,
                         rows=len(data), limits='Analytical reference only; no solver, interface or antenna acceptance.'), indent=2))


if __name__ == '__main__':
    main()
