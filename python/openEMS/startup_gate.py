# SPDX-License-Identifier: GPL-3.0-or-later
"""Controller decisions for an opt-in native startup gate; no solver launch."""
import math
import os
from pathlib import Path
import re
import tempfile
import xml.etree.ElementTree as ET


def read_startup_gate(sim_path):
    """Read a published paused-run witness, not a convergence certificate."""
    root = ET.parse(Path(sim_path) / 'startup-gate.xml').getroot()
    token = root.get('id', '')
    step = int(root.get('timestep', '0'))
    dt = float(root.get('timestepSeconds', 'nan'))
    if (root.tag != 'StartupGate' or not re.fullmatch(r'[A-Za-z0-9_-]{1,64}', token)
            or root.get('converged') != 'false' or step <= 0
            or not math.isfinite(dt) or dt <= 0):
        raise ValueError('invalid startup witness')
    return dict(id=token, timestep=step, timestep_seconds=dt, converged=False)


def decide_startup_gate(sim_path, gate_id, accept):
    """Publish accept/reject after external checks, bound to the witness ID.

    Use one controller. The filesystem must support hard links. This function
    does not establish that the solver is alive or that its physics is valid.
    """
    if type(accept) is not bool:
        raise ValueError('accept must be a boolean')
    witness = read_startup_gate(sim_path)
    if gate_id != witness['id']:
        raise ValueError('startup gate ID mismatch')
    base = Path(sim_path)
    if any((base / name).exists() for name in ['STARTUP_CONTINUE', 'STARTUP_REJECT']):
        raise FileExistsError('a startup decision already exists')
    fd, temporary = tempfile.mkstemp(prefix='.startup-decision-', dir=str(base))
    try:
        with os.fdopen(fd, 'w') as stream:
            stream.write(gate_id + '\n')
        os.link(temporary, str(base / ('STARTUP_CONTINUE' if accept else 'STARTUP_REJECT')))
    finally:
        os.unlink(temporary)
    return witness
