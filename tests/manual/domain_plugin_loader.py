# SPDX-License-Identifier: GPL-3.0-or-later
"""Reject a Metal test plugin that loads a second openEMS dylib.

Run with rebuilt Python bindings on PYTHONPATH and the explicit test plugin:
  python domain_plugin_loader.py --plugin <build>/test-plugins/<plugin>.so
No Metal device is needed; this checks the loader, not numerical execution.
"""
import argparse
import ctypes
import json
from pathlib import Path
import sys


def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument('--plugin', type=Path, required=True)
    args = parser.parse_args()
    if sys.platform != 'darwin':
        parser.error('this loader regression is specific to macOS')
    import openEMS  # load the library used by the selected Python bindings
    del openEMS
    dyld = ctypes.CDLL(None)
    dyld._dyld_image_count.restype = ctypes.c_uint32
    dyld._dyld_get_image_name.argtypes = [ctypes.c_uint32]
    dyld._dyld_get_image_name.restype = ctypes.c_char_p

    def images():
        names = [dyld._dyld_get_image_name(i).decode()
                 for i in range(dyld._dyld_image_count())]
        return sorted(name for name in names if '/libopenEMS.' in name)

    before = images()
    assert len(before) == 1, before
    plugin = ctypes.CDLL(str(args.plugin.resolve(strict=True)))
    after = images()
    assert after == before, json.dumps(dict(before=before, after=after), indent=2)
    assert plugin.openEMS_gpu_plugin, 'missing plugin entry point'
    print(json.dumps(dict(status='PASS', solver_libraries=after), indent=2))


if __name__ == '__main__':
    main()
