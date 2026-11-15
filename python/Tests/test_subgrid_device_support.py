# SPDX-License-Identifier: GPL-3.0-or-later
"""Host-only backend-log guards, not a physical GPU qualification."""
from pathlib import Path
import tempfile
import unittest
from _subgrid_device_support import run_native_device


class FakeSolver:
    def __init__(self, label, result=0, plugin='hip-nvidia'):
        self.label, self.result, self.plugin = label, result, plugin

    def Run(self, path, **options):
        print('  0: fake test GPU (' + self.plugin + ', 1.0 GiB)', flush=True)
        print('CartesianSubgrid device: ' + self.label, flush=True)
        return self.result


class NativeBackendGuards(unittest.TestCase):
    def test_no_metal_substitution_for_cuda(self):
        with tempfile.TemporaryDirectory() as path:
            with self.assertRaises(AssertionError):
                run_native_device(FakeSolver('Metal (test)'), Path(path), 'gpu', 'cuda')
            self.assertIn('Metal', (Path(path) / 'solver.log').read_text())

    def test_correct_reference_and_failed_setup(self):
        with tempfile.TemporaryDirectory() as path:
            run_native_device(FakeSolver('reference (CPU)'), Path(path) / 'reference', 'gpu-reference')
            with self.assertRaises(AssertionError):
                run_native_device(FakeSolver('HIP (test)', 4), Path(path) / 'failure', 'gpu', 'cuda')

    def test_cuda_requires_nvidia_plugin_identity(self):
        with tempfile.TemporaryDirectory() as path:
            with self.assertRaises(AssertionError):
                run_native_device(FakeSolver('HIP (test)', plugin='hip-amd'), Path(path) / 'amd', 'gpu', 'cuda')
            run_native_device(FakeSolver('HIP (test)'), Path(path) / 'nvidia', 'gpu', 'cuda')

    def test_prior_evidence_is_preserved(self):
        with tempfile.TemporaryDirectory() as path:
            run_native_device(FakeSolver('reference (CPU)'), Path(path), 'gpu-reference')
            original = (Path(path) / 'solver.log').read_bytes()
            with self.assertRaises(FileExistsError):
                run_native_device(FakeSolver('reference (CPU)'), Path(path), 'gpu-reference')
            self.assertEqual(original, (Path(path) / 'solver.log').read_bytes())

    def test_physical_selector_is_required(self):
        with tempfile.TemporaryDirectory() as path:
            with self.assertRaises(ValueError):
                run_native_device(FakeSolver('Metal (test)'), Path(path), 'gpu')


if __name__ == '__main__':
    unittest.main()
