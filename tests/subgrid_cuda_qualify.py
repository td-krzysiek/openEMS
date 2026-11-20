#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
"""Prepare/run scoped subgrid qualification on one explicitly selected CUDA GPU.

HIP is the shared backend, compiled by nvcc through OPENEMS_HIP_AS_CUDA.
This is orchestration, not another backend, scheduler or multi-GPU test.
--plan performs no hardware calls or simulation and creates no output files.
"""
import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def cache_values(path):
    result = {}
    for line in Path(path).read_text().splitlines():
        if not line or line.startswith(('//', '#')) or '=' not in line or ':' not in line.split('=', 1)[0]:
            continue
        name, value = line.split('=', 1)
        result[name.split(':', 1)[0]] = value
    return result


def enabled(value):
    return str(value).upper() in ('1', 'ON', 'YES', 'TRUE', 'Y')


def check_cache(cache):
    if not enabled(cache.get('ENABLE_HIP')):
        raise ValueError('ENABLE_HIP is OFF: configuration may have silently disabled the requested GPU backend')
    nvidia = enabled(cache.get('OPENEMS_HIP_AS_CUDA')) or cache.get('CMAKE_HIP_PLATFORM') == 'nvidia'
    if not nvidia:
        raise ValueError('NVIDIA/CUDA build required; an AMD HIP build is a separate qualification')
    if enabled(cache.get('ENABLE_METAL')) or enabled(cache.get('OPENEMS_GPU_PLUGINS')):
        raise ValueError('This scoped CUDA harness requires Metal OFF and the built-in HIP backend (plugins OFF)')
    if not enabled(cache.get('OPENEMS_BUILD_TESTS')):
        raise ValueError('OPENEMS_BUILD_TESTS=ON required')


def commands(source, build, python, include_pml=False):
    cpp = build / 'tests'
    result = [
        ('ctest', ['ctest', '--test-dir', str(build), '--output-on-failure', '-j', '2']),
        ('cuda-fields', [str(cpp / 'test_cartesian_subgrid_device'), '--cuda']),
        ('cuda-ports', [str(cpp / 'test_cartesian_subgrid_device_ports'), '--cuda']),
        ('cuda-dc', [str(cpp / 'test_cartesian_subgrid_device_ports'), '--cuda', '--dc-stress']),
        ('python-units', [python, str(source / 'python/Tests/run_testsuite.py'), '--only-unittests']),
        ('native-fields', [python, str(source / 'python/Tests/Cartesian_Subgrid_Device.py')]),
        ('native-ports', [python, str(source / 'python/Tests/Cartesian_Subgrid_Device_Ports.py')]),
    ]
    if include_pml:
        result += [
            ('cuda-pml', [str(cpp / 'test_cartesian_subgrid_device_pml'), '--cuda']),
            ('native-pml', [python, str(source / 'python/Tests/Cartesian_Subgrid_Device_PML.py')]),
            ('native-radiation', [python, str(source / 'python/Tests/Cartesian_Subgrid_Device_Radiation.py')]),
        ]
    result += [('ordinary-engines', [python, str(source / 'python/Tests/run_testsuite.py'), '--keep',
                                    'Engine_Compare', 'Cylindrical_Cavity', 'FieldProbes'])]
    return result


LOADER = r'''
import hashlib,json
from pathlib import Path
import openEMS,CSXCAD
paths=set()
for line in Path('/proc/self/maps').read_text().splitlines():
    fields=line.split()
    if len(fields)>5 and fields[-1].startswith('/'):
        p=Path(fields[-1])
        if p.is_file() and any(name in p.name for name in ['libopenEMS','libCSXCAD','libnf2ff','libcudart','libcuda.so']):
            paths.add(p.resolve())
print(json.dumps({str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},indent=2))
'''


def source_identity(source):
    names = subprocess.check_output(['git', '-C', str(source), 'ls-files', '-c', '-o',
                                     '--exclude-standard', '-z']).decode().split('\0')
    suffixes = {'.cpp', '.h', '.hpp', '.hip', '.cu', '.cuh', '.metal', '.mm', '.py', '.cmake', '.in'}
    files = {name: digest(source / name) for name in names if name and (source / name).is_file()
             and ((source / name).suffix in suffixes or name.endswith('CMakeLists.txt'))}
    return {'commit': subprocess.check_output(['git', '-C', str(source), 'rev-parse', 'HEAD'],
                                              text=True).strip(), 'files': files}


def run(args):
    source, build = args.source.resolve(), args.build.resolve()
    jobs = commands(source, build, args.python, args.include_pml)
    if args.plan:
        print(json.dumps({'status': 'PLAN_ONLY_NO_EXECUTION', 'commands': jobs,
                          'scope': 'single NVIDIA device through existing HIP backend; PML opt-in',
                          'requirements': ['qualified source/build', 'built-in NVIDIA backend', 'explicit device UUID',
                                           'CUDA driver/GPU', 'separate CPU and physical regression evidence']}, indent=2))
        return
    if not sys.platform.startswith('linux'):
        raise ValueError('The physical CUDA harness is Linux-only; use --plan on this host')
    if not args.device or not args.device.startswith('GPU-'):
        raise ValueError('Select one explicit NVIDIA UUID with --device; no automatic GPU admission')
    if args.output is None or args.output.exists():
        raise ValueError('Use a new --output directory; prior qualification must remain immutable')
    cache_path = build / 'CMakeCache.txt'
    cache = cache_values(cache_path)
    check_cache(cache)
    prefix = Path(args.prefix or cache['CMAKE_INSTALL_PREFIX']).resolve()
    if not (prefix / 'lib').is_dir():
        raise ValueError('Install the isolated candidate before qualification')
    for label, command in jobs:
        if label.startswith('cuda-') and not Path(command[0]).is_file():
            raise ValueError('Required physical test is missing: ' + command[0])
        if label.startswith('native-') and not Path(command[1]).is_file():
            raise ValueError('Required native test is missing: ' + command[1])
    output = args.output.resolve()
    output.mkdir(parents=True)
    scratch = output / 'scratch'
    scratch.mkdir()
    env = os.environ.copy()
    env['CUDA_VISIBLE_DEVICES'] = args.device
    env['OPENEMS_SUBGRID_DEVICE'] = 'cuda'
    env['LD_LIBRARY_PATH'] = str(prefix / 'lib') + ':' + env.get('LD_LIBRARY_PATH', '')
    env['TMPDIR'] = str(scratch)
    env['MPLBACKEND'] = 'Agg'
    # Log a narrow hardware/build inventory, never the full process environment.
    inventory = subprocess.check_output(['nvidia-smi', '-i', args.device,
        '--query-gpu=uuid,name,driver_version,memory.total', '--format=csv,noheader'], text=True, timeout=20)
    (output / 'nvidia-device.txt').write_text(inventory)
    loaded = json.loads(subprocess.check_output([args.python, '-c', LOADER], env=env, text=True, timeout=30))
    if not any('libopenEMS' in name and Path(name).is_relative_to(prefix) for name in loaded):
        raise ValueError('Python loaded a different core; inspect loader evidence')
    (output / 'loaded-before.json').write_text(json.dumps(loaded, indent=2) + '\n')
    source_before = source_identity(source)
    programs = {command[0]: digest(command[0]) for label, command in jobs if label.startswith('cuda-')}
    selected = {key: cache.get(key) for key in ['ENABLE_HIP', 'OPENEMS_HIP_AS_CUDA', 'CMAKE_HIP_PLATFORM',
        'CMAKE_CUDA_COMPILER', 'CMAKE_CUDA_ARCHITECTURES', 'CMAKE_CUDA_FLAGS', 'OPENEMS_GPU_BUILD_ID',
        'OPENEMS_GPU_PLUGINS', 'ENABLE_METAL', 'CMAKE_INSTALL_PREFIX']}
    manifest = {'utc': datetime.now(timezone.utc).isoformat(), 'status': 'RUNNING', 'device': args.device,
                'cacheSha256': digest(cache_path), 'cache': selected, 'source': source_before,
                'physicalTestBinaries': programs, 'commands': []}
    (output / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
    compiler = selected.get('CMAKE_CUDA_COMPILER')
    if compiler and Path(compiler).is_file():
        (output / 'cuda-compiler.txt').write_text(subprocess.check_output([compiler, '--version'], text=True, timeout=20))
    for label, command in jobs:
        start = datetime.now(timezone.utc)
        try:
            with (output / (label + '.log')).open('w') as log:
                result = subprocess.run(command, cwd=output, env=env, stdout=log, stderr=subprocess.STDOUT, timeout=3600)
        except subprocess.TimeoutExpired:
            manifest['status'] = 'TIMED_OUT'
            manifest['commands'].append({'label': label, 'command': command, 'exitCode': None,
                                         'seconds': (datetime.now(timezone.utc) - start).total_seconds()})
            (output / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
            raise RuntimeError(label + ' timed out; preserve the failed attempt')
        manifest['commands'].append({'label': label, 'command': command, 'exitCode': result.returncode,
                                     'seconds': (datetime.now(timezone.utc) - start).total_seconds()})
        manifest['status'] = 'FAILED' if result.returncode else 'RUNNING'
        (output / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
        if result.returncode:
            raise RuntimeError(label + ' failed; preserve the log and do not count as qualified')
    after = json.loads(subprocess.check_output([args.python, '-c', LOADER], env=env, text=True, timeout=30))
    if (after != loaded or digest(cache_path) != manifest['cacheSha256'] or
            source_identity(source) != source_before or any(digest(name) != sha for name, sha in programs.items())):
        raise RuntimeError('Build/runtime identity changed during qualification')
    manifest['status'] = 'PASS_SCOPED_CUDA_FIELDS_PORTS' + ('_AND_PML_DIPOLE' if args.include_pml else '')
    manifest['limits'] = ['One physical NVIDIA device; not multi-device or representative PCB antenna qualification.',
                          'No new qualification of losses/conductors or untested geometry/source regimes.']
    (output / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
    print(manifest['status'])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument('--build', type=Path, required=True)
    parser.add_argument('--prefix', type=Path)
    parser.add_argument('--python', default=sys.executable)
    parser.add_argument('--device')
    parser.add_argument('--output', type=Path)
    parser.add_argument('--include-pml', action='store_true')
    parser.add_argument('--plan', action='store_true')
    args = parser.parse_args()
    try:
        run(args)
    except (ValueError, RuntimeError, subprocess.SubprocessError, OSError) as error:
        parser.exit(1, str(error) + '\n')


if __name__ == '__main__':
    main()
