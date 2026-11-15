# SPDX-License-Identifier: GPL-3.0-or-later
"""Host-only harness checks; these do not qualify any physical CUDA kernel."""
import importlib.util
from pathlib import Path
import tempfile
import unittest

path = Path(__file__).with_name('subgrid_cuda_qualify.py')
spec = importlib.util.spec_from_file_location('subgrid_cuda_qualify', path)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


class CUDAHarness(unittest.TestCase):
    def setUp(self):
        self.valid = dict(ENABLE_HIP='ON', OPENEMS_HIP_AS_CUDA='ON', ENABLE_METAL='OFF',
                          OPENEMS_GPU_PLUGINS='OFF', OPENEMS_BUILD_TESTS='ON')

    def test_scoped_nvidia_configuration(self):
        module.check_cache(self.valid)
        for key, value in [('ENABLE_HIP', 'OFF'), ('OPENEMS_HIP_AS_CUDA', 'OFF'),
                           ('ENABLE_METAL', 'ON'), ('OPENEMS_GPU_PLUGINS', 'ON'),
                           ('OPENEMS_BUILD_TESTS', 'OFF')]:
            invalid = dict(self.valid, **{key: value})
            with self.assertRaises(ValueError, msg=key):
                module.check_cache(invalid)

    def test_physical_commands_and_explicit_pml_scope(self):
        basic = dict(module.commands(Path('/src'), Path('/build'), '/python'))
        self.assertIn('--cuda', basic['cuda-fields'])
        self.assertIn('--dc-stress', basic['cuda-dc'])
        self.assertNotIn('cuda-pml', basic)
        extended = dict(module.commands(Path('/src'), Path('/build'), '/python', True))
        self.assertIn('--cuda', extended['cuda-pml'])
        self.assertIn('native-pml', extended)
        self.assertIn('native-radiation', extended)

    def test_cache_parser(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'CMakeCache.txt'
            path.write_text('// ignored\n# comment\nENABLE_HIP:BOOL=ON\nCMAKE_CUDA_FLAGS:STRING=--fmad=false\n')
            self.assertEqual(module.cache_values(path),
                             dict(ENABLE_HIP='ON', CMAKE_CUDA_FLAGS='--fmad=false'))


if __name__ == '__main__':
    unittest.main()
