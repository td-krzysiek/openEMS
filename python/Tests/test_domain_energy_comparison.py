# SPDX-License-Identifier: GPL-3.0-or-later
"""Cross-device comparison accepts roundoff near nulls, rejects physical errors."""
import importlib.util
from pathlib import Path
import unittest

import numpy as np

source = Path(__file__).resolve().parents[2] / 'tests/manual/domain_energy.py'
spec = importlib.util.spec_from_file_location('domain_energy_comparison', source)
comparison = importlib.util.module_from_spec(spec)
spec.loader.exec_module(comparison)


class DomainEnergyComparison(unittest.TestCase):
    def test_near_null_roundoff_uses_reference_peak(self):
        reference = np.array([1 + 1j, 1e-12, 0])
        candidate = reference + np.array([1e-7, 1e-8, 0])
        self.assertLess(comparison.peak_relative_error(candidate, reference), 1e-6)

    def test_significant_error_rejected(self):
        with self.assertRaises(AssertionError):
            comparison.peak_relative_error(np.array([1.001]), np.array([1.0]))

    def test_zero_reference_requires_exact_zero(self):
        self.assertEqual(comparison.peak_relative_error(np.zeros(2), np.zeros(2)), 0)
        with self.assertRaises(AssertionError):
            comparison.peak_relative_error(np.array([0, 1e-30]), np.zeros(2))

    def test_invalid_data_rejected(self):
        for candidate in [np.array([np.nan]), np.array([np.inf]), np.ones(2)]:
            with self.subTest(candidate=candidate):
                with self.assertRaises(AssertionError):
                    comparison.peak_relative_error(candidate, np.ones(1))


if __name__ == '__main__':
    unittest.main()
