# SPDX-License-Identifier: GPL-3.0-or-later
"""Atomic run-control requests; no native solver or GPU required."""
import importlib.util
from pathlib import Path
import tempfile
import unittest

spec = importlib.util.spec_from_file_location('run_control', Path(__file__).parents[1] / 'openEMS/run_control.py')
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


class RunControlTest(unittest.TestCase):
    def test_atomic_pause_refuses_pending_request(self):
        with tempfile.TemporaryDirectory() as folder:
            module.control_run(folder, 'pause', 'run-manual')
            self.assertEqual((Path(folder) / 'RUNTIME_PAUSE').read_text(), 'run-manual\n')
            with self.assertRaises(FileExistsError):
                module.control_run(folder, 'pause', 'run-other')
            self.assertEqual(len(list(Path(folder).iterdir())), 1)

    def test_decision_matches_paused_id(self):
        with tempfile.TemporaryDirectory() as folder:
            (Path(folder) / 'runtime-gate.xml').write_text('<RuntimeGate id="run-manual" state="paused"/>')
            with self.assertRaises(ValueError):
                module.control_run(folder, 'continue', 'another')
            module.control_run(folder, 'continue', 'run-manual')
            self.assertEqual((Path(folder) / 'RUNTIME_CONTINUE').read_text(), 'run-manual\n')

    def test_invalid_inputs_do_not_write(self):
        with tempfile.TemporaryDirectory() as folder:
            for action, token in [('exec', 'run'), ('pause', '../run'), ('stop', 'x' * 65)]:
                with self.assertRaises(ValueError):
                    module.control_run(folder, action, token)
            self.assertEqual(list(Path(folder).iterdir()), [])


if __name__ == '__main__':
    unittest.main()
