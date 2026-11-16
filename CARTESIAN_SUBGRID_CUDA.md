# Preparing later Cartesian-subgrid tests on NVIDIA/CUDA

This is a build/reproduction note. The canonical status is
[CARTESIAN_SUBGRID_CHECKLIST.md](CARTESIAN_SUBGRID_CHECKLIST.md), not this file.
Physical Metal controls and the separately pinned single-RTX3080 lane below
have scoped evidence. Merely preparing CUDA code/host harness tests does not
qualify another source revision, NVIDIA device or multi-device subgridding.

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


## Qualified single RTX3080 lane

Numerical source: `02b94efcd1f5b2546569f45e60a246285ea3041b`, with only the
explicit `<cstdint>` test-header correction later committed as `a846ebf`.
The collected source inventory independently matches those public git blobs;
no numerical engine difference is allowed. Matching CSXCAD headers/library and
fresh Python bindings use `cc37cb4ee0ef20c6794acc3361b5f56953eb5189`.

Existing HIP is built into the core using NVIDIA nvcc 12.6.85, HIP NVIDIA headers,
Metal/plugins OFF, SM86/89 compiled and `--fmad=false`. Physical execution uses
one RTX3080 (SM86), UUID `GPU-d9f1434a-b4b4-c89e-658f-83f93d632df0`.
Compiling SM89 does not qualify an RTX4090. The tests verify the actual HIP
backend and loaded libraries rather than accepting CPU/reference fallback.

All eleven original `--include-pml` harness gates pass across four immutable
attempts: CTest 21/21, physical fields/ports/DC, Python 14/14, native fields/ports,
core/native PML, public vacuum PEC dipole/host NF2FF and three ordinary-engine
groups. The radiation gate includes identity/local meshes, complex near/far
fields, two closed surfaces, accepted/radiated/sphere power, source contacts,
TD/FD and wrong-half-step/missing-face negatives. Fields, port signs and native
input checks retain their original tolerances.

| Quantity | Maximum measured error | Gate |
| --- | ---: | ---: |
| Physical full fields vs double CPU | 8.975e-6 | 2e-4 |
| Physical fields vs independent FP32 oracle | 0 | 2e-6 |
| Strong-DC fields | 1.761e-5 | 2e-4 |
| Native port complex impedance | 1.004e-5 ohm | 0.01 ohm |
| Dipole complex impedance | 4.927e-6 ohm | 0.01 ohm |
| Dipole normalized field/power comparisons | 4.914e-7 | 2e-4 |

Attempts r01/r02 hit external 20-minute deadlines, without OOM. Attempt r03
passed ordinary engines but the radiation child exceeded our 1100-second
concurrent orchestration limit. These remain incomplete/timeout records, not
relabeled completed runs. The final attempt verified their terminal proofs and
reran only the unchanged missing radiation gate, alone with 4 CPU / 3 GiB host RAM
and the original harness 3600-second per-gate limit (3700s outer bound). It
completed in 1208.93 seconds. All eleven gate labels and final source, cache,
test-binary, orchestrator and actual loaded-library identities pass. No numerical
failure was retried as though it were a timeout.

The shared-card test was explicitly authorized alongside an active antenna;
its resource sampler could stop only its own container. No active model,
production runtime, lease or unrelated process was changed. These concurrent
timings are not an uncontended CUDA throughput benchmark.

Loaded hashes:

- Core: `6f711fd7afd10ea6c758cfee09ae4f5f3fefaa2aaa6e8690e4bd1ba87fdcac9a`.
- NF2FF: `22cf200dde5d2ca350d3d045c424bcd256e3be74df25146806989cd0adebaad9`.
- CSXCAD: `7482cb5a7a86e447ee4b960224902e8e0ec819e4747a205af48e225f6266275b`.

Contributor evidence is under
`experiments/openems-cartesian-local-subgrid-20261112/cuda/`: final
`qualification.json`, four manifests/logs, source/runtime inventories and a
4.38 MB compact archive; raw TD faces remain in the isolated remote test root.
Earlier build failures are retained: mismatched old CSXCAD headers, then the
GCC test include. Neither is presented as an upstream numerical defect.

This qualifies only the stated pinned single-device lossless/PEC field, source,
one-edge port, PML and dipole scope. The later native polygon and larger-interface
increments need fresh physical CUDA acceptance. Plugin mode, RTX4090 hardware,
lossy/surface/ADE conductors, representative PCB accuracy and multiple devices
remain open. Existing production multi-GPU execution does not qualify this
compact subgrid path.

## Scoped logical domains on one RTX3080

The [domain increment](CARTESIAN_SUBGRID_DOMAINS.md), configure-time parent
`1a8f48f` plus its exact recorded source diff, uses the existing built-in HIP
backend compiled by nvcc 12.6.85 for SM86/SM89, with `--fmad=false`. Only physical
RTX3080 SM86 is tested. The same numerical source passes reference/Metal locally.

Run `test_cartesian_subgrid_device_domains --cuda` from the matching build. It
rejects fallback, exercises 1–4 logical domains and mixed CPU-reference/CUDA
workers, and retains the original field gates, phase/epoch and invalid-input
controls. All eight 1500-step cases are identical to one FP32 device, with
maximum double error 1.62522e-6. CTest 23/23 and existing single-device fields,
ports and DC controls also pass; maximum physical field error is 8.975e-6 and
DC error 1.761e-5 (gate 2e-4). These port controls are not domain-port support.

The isolated container exposes one explicit GPU UUID, 4 host CPUs and 4 GiB RAM.
The first attempt reached our external 300 s CPU CTest deadline after 20 passes;
it remains incomplete evidence. The second uses identical source/runtime/gates,
runs physical controls first and allows 900 s for CTest. It completes in 724 s,
including 647 s of CPU CTest; no concurrent-throughput claim is made. Final
container exit is 0/noOOM and all source/runtime identities remain unchanged.

CUDA core SHA256:
`de9994d3e5a09f07de74de1effed5cc9c6371d2d7ee28bd0204d8ccdeab56c0a`.
Compact 54-file evidence and the exact source/cache/test identities are retained
in the project study's `domains/cuda-qualification.json` and
`output/openems-subgrid-domains-cuda-20261115/completed-r02/`. This does not
qualify physical multi-GPU, native CAD/large-interface CUDA, plugin mode,
distributed ports/PML/DFT, losses or representative PCB antenna workflows.

## Coupled-domain-port CUDA lane

The [port-domain increment](CARTESIAN_SUBGRID_DOMAINS.md#coupled-resistors-and-source-ledgers),
configure-time parent `15daf56` plus its recorded source diff, adds complete small
port-solve ownership and per-worker source/loss accounting. Six three-worker
cases pass on one RTX3080, including mixed CPU-reference/CUDA, source-only
workers, warm accumulated ledgers and strong DC. Fields are identical to one
FP32 device; maximum double error 1.94051e-5 and balance 1.70516e-6 (gates 2e-4).

All six campaign gates pass: retained bare domains, coupled port domains,
fields, ports, DC and CTest 24/24. The first bounded attempt completes and exits 0
without OOM; no numerical gate or model changed during the campaign. CTest
uses 587 s in this concurrent run; no throughput claim. Final source/cache/binary
checks pass. The 48-file compact archive and independent collector are in the
project study's `domain-ports/` and
`output/openems-subgrid-domain-ports-cuda-20261116/completed-r01/`.

Core SHA256: `f9af217240707dba76822509402fd835c5a5a47cbbf10e4eef121d79096ef0de`.
Built-in HIP/nvcc 12.6.85, SM86/SM89 compiled, only RTX3080 SM86 physically tested.
Physical multi-GPU, native domain RF/processing, PML/loss histories, plugin mode
and representative PCB/conductor accuracy remain unqualified.
