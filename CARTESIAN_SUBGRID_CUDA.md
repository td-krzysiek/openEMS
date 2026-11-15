# Preparing later Cartesian-subgrid tests on NVIDIA/CUDA

This is a build/reproduction note. The canonical status is
[CARTESIAN_SUBGRID_CHECKLIST.md](CARTESIAN_SUBGRID_CHECKLIST.md), not this file.
The current physical qualification is local Metal. Preparing CUDA code and
host-only harness tests does **not** qualify NVIDIA hardware or multi-device
subgridding.

## Same HIP backend, NVIDIA compiler/runtime

Use the existing HIP backend and shared portable point operations. With
`OPENEMS_HIP_AS_CUDA=ON`, CMake compiles the existing `.hip` sources as CUDA
using `nvcc`. HIP's NVIDIA headers map the API to CUDA; an AMD ROCm runtime
or a separate new backend is not required. CUDA toolkit/driver and the HIP
headers are still build dependencies; see the pinned master's CMake discovery.
Preserve `--fmad=false`/safe arithmetic except for the deliberately explicit
FMA operations in the compensated subgrid code. Do not enable global fast math
to make this configuration compile or appear faster.

Use an isolated build/install prefix and a source checkout without the local
Mac-only `localConfig.cmake` overrides. Do not change or replace a running or
qualified runtime. Example configuration for RTX3080 (`sm_86`) and RTX4090
(`sm_89`), with existing dependency prefixes supplied as appropriate:

```sh
cmake -S . -B build-subgrid-cuda \
  -DENABLE_HIP=ON -DOPENEMS_HIP_AS_CUDA=ON -DENABLE_METAL=OFF \
  -DOPENEMS_GPU_PLUGINS=OFF -DOPENEMS_BUILD_TESTS=ON \
  '-DCMAKE_CUDA_ARCHITECTURES=86;89' \
  -DCMAKE_INSTALL_PREFIX=/your/isolated/subgrid-cuda-prefix \
  -DCSXCAD_ROOT_DIR=/your/dependency-prefix \
  -DHIP_INCLUDE_DIR=/your/HIP-headers/include
cmake --build build-subgrid-cuda -j 8
cmake --install build-subgrid-cuda
```

Check `CMakeCache.txt` and the configure log: missing HIP headers/CUDA tools may
cause the upstream configuration to disable the requested backend. A successful
CPU-only build must never be described as CUDA qualification. This first harness
intentionally tests the built-in HIP backend (`OPENEMS_GPU_PLUGINS=OFF`);
plugin-mode coverage is a separate lane. The C++ `--cuda` selector exists only
in a test built for that NVIDIA configuration, selects device0 and verifies an
actual HIP backend. `--metal` still verifies Metal. Neither request falls back
to the CPU reference.

## Plan first, then physical evidence

The runner uses only the Python standard library for orchestration. Native
solver checks still need the matching openEMS/CSXCAD Python bindings and their
existing NumPy/HDF5 dependencies. Set `PYTHONPATH`/dependency library paths for
that isolated build, and select an admitted physical GPU explicitly by UUID.
The runner is not a scheduler and does not preempt another process.

```sh
python3 tests/subgrid_cuda_qualify.py --build build-subgrid-cuda --plan
python3 tests/subgrid_cuda_qualify.py --build build-subgrid-cuda \
  --device GPU-YOUR-EXACT-UUID --output /your/new/cuda-evidence-directory
```

`--plan` performs no hardware calls, simulations or output-directory writes;
it can run on macOS. Physical execution is Linux-only, uses a fresh evidence
directory and one `CUDA_VISIBLE_DEVICES` UUID. It requires the expected CUDA
configuration, records the selected NVIDIA/driver information, compiler,
actual loaded core/dependencies, source/test-binary/cache identities, command
results and logs. Native test logs independently reject Metal/reference backend
substitution when CUDA was requested.
CTest/default reference, physical fields/ports/DC, native scripting and ordinary
engine controls all run. A missing GPU, wrong build, failed gate or changed
runtime fails the campaign rather than skipping it.

The default scope is closed-boundary fields and ports. `--include-pml` also
requires the separate integrated device-PML test binary/native script and
public PEC dipole port/near-field/NF2FF controls; it
must fail if these are not yet built/present. It must not imply that PML or an
antenna was qualified merely because the basic field/port lane passed.
Compare immutable results with the same-input CPU/Metal controls, including
complex fields and S parameters, energy/work, and PML/radiation where requested.
Record device memory and completed-batch times separately from setup/DFT.
Only after those numerical gates pass should a representative antenna or
multi-GPU experiment be planned. Existing upstream CUDA/HIP support alone
does not complete these subgrid gates.
