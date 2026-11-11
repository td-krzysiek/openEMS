# SPDX-License-Identifier: GPL-3.0-or-later
"""Build and run the disjoint-region regression using a variant's include paths."""
import argparse
import os
from pathlib import Path
import shlex
import subprocess

parser = argparse.ArgumentParser()
parser.add_argument('--build', type=Path, required=True)
parser.add_argument('--prefix', type=Path, required=True)
parser.add_argument('--test', choices=['test_upml_disjoint.cpp', 'test_material_dispatch.cpp', 'test_iteration_barrier.cpp'],
                    default='test_upml_disjoint.cpp')
args = parser.parse_args()
flags = (args.build / 'CMakeFiles/openEMS.dir/flags.make').read_text().splitlines()
includes = shlex.split(next(line.split('=', 1)[1] for line in flags if line.startswith('CXX_INCLUDES =')))
binary = args.build / Path(args.test).stem
command = ['c++', '-std=c++11', '-O2', '-pthread', *includes,
           str(Path(__file__).with_name(args.test)), '-o', str(binary),
           '-L' + str(args.build), '-L' + str(args.prefix / 'lib'), '-lopenEMS', '-lCSXCAD',
           '-Wl,-rpath-link,' + str(args.prefix / 'lib')]
subprocess.run(command, check=True, timeout=120)
environment = dict(os.environ, LD_LIBRARY_PATH=f'{args.build}:{args.prefix / "lib"}')
subprocess.run([str(binary)], env=environment, check=True, timeout=30)
