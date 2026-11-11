# SPDX-License-Identifier: GPL-3.0-or-later
"""Controller boundary checks; native continuation has a separate manual test."""
import importlib.util
from pathlib import Path
import tempfile
import unittest

source = Path(__file__).resolve().parents[1] / 'openEMS/startup_gate.py'
spec = importlib.util.spec_from_file_location('startup_controller', source)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


class StartupDecision(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.path = Path(self.directory.name)
        (self.path / 'startup-gate.xml').write_text(
            '<StartupGate id="current" timestep="4100" timestepSeconds="1e-14" converged="false"/>')

    def test_current_decision_is_exclusive(self):
        result = module.decide_startup_gate(self.path, 'current', True)
        self.assertEqual(result['timestep'], 4100)
        self.assertEqual((self.path / 'STARTUP_CONTINUE').read_text(), 'current\n')
        with self.assertRaises(FileExistsError):
            module.decide_startup_gate(self.path, 'current', False)
        self.assertFalse((self.path / 'STARTUP_REJECT').exists())

    def test_foreign_id_cannot_release_gate(self):
        with self.assertRaises(ValueError):
            module.decide_startup_gate(self.path, 'old', True)
        self.assertFalse((self.path / 'STARTUP_CONTINUE').exists())

    def test_invalid_witness_does_not_publish(self):
        (self.path / 'startup-gate.xml').write_text(
            '<StartupGate id="current" timestep="4100" timestepSeconds="nan" converged="false"/>')
        with self.assertRaises(ValueError):
            module.decide_startup_gate(self.path, 'current', True)
        self.assertFalse((self.path / 'STARTUP_CONTINUE').exists())

    def test_reject_remains_reject(self):
        module.decide_startup_gate(self.path, 'current', False)
        self.assertEqual((self.path / 'STARTUP_REJECT').read_text(), 'current\n')


if __name__ == '__main__':
    unittest.main()
