# SPDX-License-Identifier: GPL-3.0-or-later
"""Generic in-memory pause control for an opt-in openEMS process.

Use --runtime-control --decay-gate-id=RUN_TOKEN. Requests are acted on at a
processing boundary. These files do not restart a process after it exits.
"""
import os
from pathlib import Path
import re
import tempfile
import xml.etree.ElementTree as ET


def control_run(sim_path, action, request_id):
    """Atomically request ``pause``, ``continue`` or graceful ``stop``.

    A manual pause ID must start with the launch's run token followed by '-'.
    Read runtime-gate.xml to obtain the actual paused step; continue/stop must
    match its ID. The existing startup protocol remains backward compatible.
    """
    if action not in ('pause', 'continue', 'stop'):
        raise ValueError('action must be pause, continue or stop')
    if not isinstance(request_id, str) or not re.fullmatch(r'[A-Za-z0-9_-]{1,64}', request_id):
        raise ValueError('Invalid request ID')
    base = Path(sim_path)
    if action != 'pause':
        witness = ET.parse(base / 'runtime-gate.xml').getroot()
        if witness.tag != 'RuntimeGate' or witness.get('id') != request_id or witness.get('state') != 'paused':
            raise ValueError('Decision does not match the paused runtime gate')
    fd, temporary = tempfile.mkstemp(prefix='.run-control-', dir=str(base))
    try:
        with os.fdopen(fd, 'w') as stream:
            stream.write(request_id + '\n')
        os.link(temporary, base / ('RUNTIME_' + action.upper()))
    finally:
        os.unlink(temporary)
