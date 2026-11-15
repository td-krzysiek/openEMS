# SPDX-License-Identifier: GPL-3.0-or-later
import json
import math


def read_setup_profile(text):
    """Read the optional native setup record, retaining inclusive phase hierarchy."""
    records = [line[len('Setup profile: '):] for line in text.splitlines()
               if line.startswith('Setup profile: ')]
    if not records:
        return None
    if len(records) != 1:
        raise ValueError('expected one setup profile per solver invocation')
    profile = json.loads(records[0])
    if not isinstance(profile, dict) or profile.get('schema') != 'openems.setup-profile.v1':
        raise ValueError('unsupported setup profile schema')
    if profile.get('status') not in ('complete', 'incomplete'):
        raise ValueError('invalid setup profile status')
    if profile.get('scope') not in ('setup', 'operator-only'):
        raise ValueError('invalid setup profile scope')
    def valid_seconds(value):
        return type(value) in (int, float) and math.isfinite(value) and value >= 0
    if not valid_seconds(profile.get('seconds')) or not isinstance(profile.get('phases'), list):
        raise ValueError('invalid setup profile timing')
    for index, phase in enumerate(profile['phases']):
        if not isinstance(phase, dict):
            raise ValueError('invalid setup profile phase')
        parent = phase.get('parent')
        if (type(parent) is not int or not -1 <= parent < index
                or not isinstance(phase.get('name'), str)
                or not isinstance(phase.get('detail'), str)
                or not valid_seconds(phase.get('seconds'))):
            raise ValueError('invalid setup profile phase')
    return profile
