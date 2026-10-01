"""Transient tiling regressions. Run against a built/installed candidate.

Each solver runs in a subprocess with a timeout, so a worker shutdown deadlock
fails the test. The tolerance is 1e-6 of the reference array's peak plus 1e-20;
this checks numerical equivalence, not antenna convergence or performance.
"""
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import numpy as np


def fixture(destination, engine, threads, case):
    from CSXCAD import ContinuousStructure
    from openEMS import openEMS
    import h5py

    destination = Path(destination).resolve()
    destination.mkdir(parents=True)
    fdtd = openEMS(NrTS=1203, EndCriteria=0, OverSampling=4,
                   CoordSystem=1 if case == 'cylinder' else 0)
    fdtd.SetGaussExcite(8e9, 6e9)
    fdtd.SetBoundaryCond(['MUR'] * 6 if case == 'mur' else ['PML_8'] * 6)
    csx = ContinuousStructure()
    fdtd.SetCSX(csx)
    mesh = csx.GetGrid()
    mesh.SetDeltaUnit(1e-3)
    # Odd sizes, uneven Z spacing, partial XY tiles and padded SIMD rows.
    mesh.SetLines('x', np.linspace(-24, 24, 65))
    mesh.SetLines('y', np.linspace(-21, 21, 57))
    mesh.SetLines('z', np.unique(np.r_[np.linspace(-15, 18, 45), 0, 0.5, 1]))
    csx.AddMaterial('substrate', epsilon=3.66, kappa=0.002).AddBox(
        [-15, -12, 0], [15, 12, 1])
    if case == 'dispersive':
        from CSXCAD.CSProperties import CSPropLorentzMaterial
        material = CSPropLorentzMaterial(csx.GetParameterSet(), name='lorentz',
                                        epsilon=1.8, mue=1.2, order=2)
        csx.AddProperty(material)
        material.SetDispersiveMaterialProperty(0, eps_plasma=3e9, eps_pole_freq=12e9,
                                              eps_relax=50e-12, mue_plasma=1e9,
                                              mue_pole_freq=10e9, mue_relax=60e-12)
        material.SetDispersiveMaterialProperty(1, eps_plasma=2e9, eps_pole_freq=0,
                                              eps_relax=80e-12, mue_plasma=0.5e9,
                                              mue_pole_freq=0, mue_relax=80e-12)
        material.AddBox([-6, -6, 3], [6, 6, 6])
    csx.AddMetal('ground').AddBox([-15, -12, 0], [15, 12, 0], priority=10)
    csx.AddConductingSheet('patch', conductivity=5.8e7, thickness=35e-6).AddBox(
        [-4.5, -3.75, 1], [4.5, 3.75, 1], priority=10)
    ports = [fdtd.AddLumpedPort(1, 50, [-1.5, 0, 0], [-1.5, 0, 1],
                               'z', excite=1, priority=20),
             fdtd.AddLumpedPort(2, 75, [3, 0, 0], [3, 0, 1],
                               'z', excite=0, priority=20)]
    if case == 'rlc':
        csx.AddLumpedElement('load', ny=2, caps=True, R=120, L=1e-9, C=0.2e-12).AddBox(
            [0, 2.25, 0], [0, 2.25, 1], priority=20)
    for name, kind in [('electric', 10), ('magnetic', 11)]:
        csx.AddDump(name, dump_type=kind, file_type=1, frequency=[8e9]).AddBox(
            [-18, -15, -9], [18, 15, 12])
    status = fdtd.Run(str(destination), cleanup=False, numThreads=threads, engine=engine)
    if status:
        # Exercise destruction after rejected initialization, too.
        fdtd.Reset()
        raise SystemExit(10)
    frequency = np.linspace(3e9, 13e9, 101)
    arrays = {'frequency': frequency}
    for i, port in enumerate(ports, 1):
        port.CalcPort(str(destination), frequency, ref_impedance=50)
        arrays[f'incident{i}'] = port.uf_inc
        arrays[f'reflected{i}'] = port.uf_ref
    for path in sorted(destination.glob('port_*')):
        if path.is_file():
            arrays[path.name] = np.loadtxt(path, comments='%')
    for name in ['electric', 'magnetic']:
        with h5py.File(destination / (name + '.h5')) as field:
            def visit(key, obj):
                if isinstance(obj, h5py.Dataset):
                    arrays[name + '/' + key] = np.asarray(obj[()])
            field.visititems(visit)
    np.savez(destination / 'arrays.npz', **arrays)
    maps = Path('/proc/self/maps')
    if maps.exists():
        (destination / 'loaded-libraries.txt').write_text('\n'.join(
            row for row in maps.read_text().splitlines() if 'libopenEMS' in row))
    fdtd.Reset()  # explicit worker teardown, followed by destructor at exit


class TestTiling(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        requested = os.environ.get('OPENEMS_TILING_TEST_OUTPUT')
        cls.temp = None if requested else tempfile.TemporaryDirectory()
        cls.root = Path(requested or cls.temp.name).resolve()
        cls.root.mkdir(parents=True, exist_ok=True)
        cls.report = []

    @classmethod
    def tearDownClass(cls):
        (cls.root / 'comparison.json').write_text(json.dumps(cls.report, indent=2) + '\n')
        if cls.temp:
            cls.temp.cleanup()

    def run_case(self, name, case='pml', engine='tiling', threads=2, tiles=(10, 10, 5), rejection=None):
        destination = self.root / name
        env = dict(os.environ)
        for key, value in zip(['X', 'Y', 'STEPS'], tiles):
            env['OPENEMS_TILE_' + key] = str(value)
        proc = subprocess.run([sys.executable, str(Path(__file__).resolve()), '--fixture',
                               str(destination), engine, str(threads), case],
                              env=env, capture_output=True, text=True, timeout=90)
        (self.root / (name + '.log')).write_text(proc.stdout + proc.stderr)
        if rejection:
            self.assertEqual(proc.returncode, 10, proc.stdout + proc.stderr)
            self.assertIn(rejection, proc.stdout + proc.stderr)
            return
        self.assertEqual(proc.returncode, 0, proc.stdout + proc.stderr)
        if engine == 'tiling':
            self.assertIn('Experimental tiling:', proc.stdout + proc.stderr)
        with np.load(destination / 'arrays.npz') as data:
            return {key: data[key] for key in data.files}

    def compare(self, name, reference, candidate):
        self.assertEqual(set(reference), set(candidate))
        self.assertGreater(np.max(np.abs(reference['port_ut_1'][:, 1])), 1e-4)
        for key, expected in reference.items():
            actual = candidate[key]
            self.assertEqual(expected.shape, actual.shape)
            self.assertTrue(np.all(np.isfinite(expected)) and np.all(np.isfinite(actual)))
            peak = float(np.max(np.abs(expected)))
            error = float(np.max(np.abs(actual - expected)))
            passed = error <= 1e-6 * peak + 1e-20
            self.report.append(dict(case=name, array=key, peak=peak, max_error=error,
                                    identical=bool(np.array_equal(expected, actual)), passed=passed))
            self.assertTrue(passed, f'{name}/{key}: error {error}, peak {peak}')

    def test_transients(self):
        for case in ['pml', 'dispersive']:
            ref = self.run_case(case + '-ref', case, engine='multithreaded')
            # Single-step spatial tiling, temporal batches, non-divisible batches,
            # tiles exceeding the domain, and multiple worker counts.
            for threads, tiles in [(1, (10, 10, 5)), (2, (7, 11, 1)),
                                   (4, (10, 10, 5)), (4, (3, 5, 8)), (2, (100, 100, 2))]:
                name = f'{case}-{threads}-' + '-'.join(map(str, tiles))
                with self.subTest(name=name):
                    self.compare(name, ref, self.run_case(name, case, threads=threads, tiles=tiles))

    def test_rejections(self):
        for case in ['mur', 'rlc']:
            with self.subTest(case=case):
                self.run_case('reject-' + case, case, rejection='does not support the experimental tiling engine')
        self.run_case('reject-zero', tiles=(0, 10, 5), rejection='OPENEMS_TILE_X')
        self.run_case('reject-negative', tiles=(10, 10, -1), rejection='OPENEMS_TILE_STEPS')
        self.run_case('reject-cylinder', 'cylinder', rejection='does not support cylindrical')


if __name__ == '__main__':
    if len(sys.argv) > 1 and sys.argv[1] == '--fixture':
        fixture(sys.argv[2], sys.argv[3], int(sys.argv[4]), sys.argv[5])
    else:
        unittest.main()
