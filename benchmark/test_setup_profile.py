# SPDX-License-Identifier: GPL-3.0-or-later
import copy
import json
import unittest
from setup_profile import read_setup_profile


class SetupProfileTest(unittest.TestCase):
    def setUp(self):
        self.profile = dict(schema='openems.setup-profile.v1', status='complete',
                            scope='setup', seconds=2.0,
                            phases=[dict(name='operator', detail='', parent=-1, seconds=1.9),
                                    dict(name='operator.extension', detail='PML', parent=0, seconds=0.2)])

    def parse(self, value):
        return read_setup_profile('existing log\nSetup profile: ' + json.dumps(value) + '\n')

    def test_optional(self):
        self.assertIsNone(read_setup_profile('Creation time for operator: 2 s\n'))

    def test_hierarchy_and_repeated_names(self):
        self.profile['phases'].append(copy.deepcopy(self.profile['phases'][1]))
        self.assertEqual(self.parse(self.profile), self.profile)

    def test_reject_invalid(self):
        for key, value in [('schema', 'unknown'), ('status', 'wrong'), ('scope', 'wrong'),
                           ('seconds', -1), ('seconds', float('nan')), ('seconds', True)]:
            with self.subTest(key=key, value=value):
                profile = dict(self.profile, **{key: value})
                with self.assertRaises(ValueError):
                    self.parse(profile)
        for key, value in [('parent', 1), ('parent', True), ('seconds', float('inf')), ('detail', None)]:
            with self.subTest(key=key, value=value):
                profile = copy.deepcopy(self.profile)
                profile['phases'][1][key] = value
                with self.assertRaises(ValueError):
                    self.parse(profile)

    def test_malformed_or_multiple(self):
        with self.assertRaises(ValueError):
            read_setup_profile('Setup profile: {oops}\n')
        record = 'Setup profile: ' + json.dumps(self.profile) + '\n'
        with self.assertRaises(ValueError):
            read_setup_profile(record * 2)


if __name__ == '__main__':
    unittest.main()
