# Experimental domain runtime

Internal C++11 building blocks for one Cartesian FDTD problem across compute
resources. This is not an antenna-ready engine or a new `--engine` selection.
Existing CPU/GPU engines and Python/Matlab interfaces remain unchanged.

## Code style

Use explicit braces for every conditional and loop, including one-statement
bodies. Expand functions and control flow over multiple lines. The scoped
`.clang-format` configuration uses Allman braces, four-column tab indentation
and a 100-column limit (checked with clang-format 23.1.3).

Apply the same configuration to the executor and test files:

```sh
clang-format --style=file:FDTD/distributed/.clang-format -i \
  FDTD/distributed/*.{h,cpp} FDTD/hip/hip_domain_executor.{h,hip} \
  TESTSUITE/distributed/*.{h,cpp}
```

## Boundaries

| Layer | Responsibility |
|---|---|
| `PartitionPlan` | Unique native-node ownership, local/ghost maps and directional tangential halos |
| Placement | Domain-to-host/process/backend/resource mapping and aggregate memory budget |
| `RunManager` | Startup validation, lifecycle and a single logical run |
| `FDTDCoordinator` | Bootstrap, corrected voltage/current exchanges and phase-specific observations |
| `DomainExecutor` | Local updates, packing/unpacking and their completion tokens |
| `HaloTransport` | Local mailbox or optional MPI movement; no numerical interpretation |
| `CollectiveContext` | Global timestep agreement, ownership checks, reductions and failure propagation |

The same core runs with several executors in one process or one coordinator
per MPI rank. A rank may own multiple domains. No manager relays halo fields.
MPI uses nonblocking header/payload transfers, explicit progress on the MPI
initializing thread, bounded tags and separate communicators. The header
validates plan, step, phase, channel, endpoints, count and FP32 format. A failed
rank aborts the MPI run; recovery is not implemented.

Receive-ready includes destination copying and unpacking. Send-release is a
different completion. Exchange buffers are reused only after both complete.
The initial implementation waits at half-step boundaries; it does not overlap
interior work or use production full-timestep fusion. Adapters use no virtual
dispatch inside a cell-update loop.

## Supported scope

- Scalar CPU and HIP domain executors, fixed nonempty X slabs, global dt,
  static final operator coefficients and prescribed initial fields.
- Actual openEMS operator slicing through `SliceOperator`, preserving global
  averages and true boundaries. Cartesian PEC/static dielectric coefficients
  are supported. All active dynamic extensions except an empty excitation
  extension are rejected before updates.
- Pinned-host HIP staging. The portable host buffer currently adds a CPU copy
  at each staging endpoint. No P2P or GPU-aware MPI claim is made.
- Optional CPU-only MPI executable can communicate with a separate HIP
  executable without linking the CPU executable to a GPU runtime.

Not implemented: distributed UPML, opaque poles, driven/lumped ports, complex
DFT/NF2FF, production energy/stop processing, Metal executor, online balancing
or state migration. The observation hook is a phase-specific contract, not a
replacement for existing production processing. Real-antenna inputs containing
these extensions must fail; no physics is silently dropped.

Fields and coefficient arrays are local, with two halo storage planes. The
basic uncompressed estimate is six field arrays over owned+ghost storage plus
twelve owned coefficient arrays. Auxiliary-state estimates must be supplied
separately when supported. Shared physical memory pools use the same host and
resource key regardless of backend; budgets aggregate across their domains.
Driver/context, staging and construction peaks require additional preflight.
Building the global operator once does not reduce peak host setup RAM.

## Build and test

Neutral core tests need neither openEMS dependencies nor GPU/MPI:

```sh
cmake -S FDTD/distributed -B build-domain -DOPENEMS_DOMAIN_TESTS=ON
cmake --build build-domain
ctest --test-dir build-domain --output-on-failure
```

Add `-DOPENEMS_DOMAIN_MPI=ON` where MPI is already available. Tests cover one,
two and three ranks, multiple domains per rank and a wrong-epoch abort control.

For operator/HIP integration, configure the full project with its normal
dependency prefixes and these opt-in settings:

```text
ENABLE_DOMAIN_RUNTIME=ON
OPENEMS_DOMAIN_TESTS=ON
OPENEMS_DOMAIN_HIP=ON
OPENEMS_DOMAIN_MPI=ON
```

The HIP option requires the existing qualified HIP/CUDA toolchain. Build
`test_domain_hip`, `test_domain_mpi_hip`, `test_domain_mpi_cpu` and
`test_domain_operator`. Expose exactly two idle devices for `test_domain_hip`;
the test driver uses their runtime ordinals and prints PCI mappings.

`TESTSUITE/distributed/make_operator_fixtures.py <directory>` writes small
PEC/dielectric XML controls using installed Python bindings.
`test_domain_operator fixture.xml snapshot.txt` checks four-domain CPU fields
against the basic engine and exports a bounded coefficient test fixture.
Pass it to `test_domain_hip --fixture snapshot.txt` or to the MPI executables.
The fixture loader is test-only; it intentionally bounds input size.

```sh
mpirun -np 2 ./test_domain_mpi_hip snapshot.txt
mpirun -np 1 ./test_domain_mpi_cpu snapshot.txt --cpu-gpu : \
       -np 1 ./test_domain_mpi_hip snapshot.txt --cpu-gpu
```

These tests use a peak-normalized six-field gate of `2e-5`; wrong/missing
halos must fail their negative controls. Same-arithmetic CPU split cases also
check exact agreement. `test_domain_hip --benchmark` times 512 steps after
32 warmup steps, with three reversed-order repetitions. It compares this
unfused runtime with itself, not the fastest production GPU engine; it cannot
establish antenna speedup. Full models/measurements are not committed here.

## Initial qualification (20 October 2026)

- C++11 core and CPU/MPI controls pass on Linux; macOS core tests also pass
  ASan/UBSan. The complete macOS CPU/Metal project builds with the option on.
- Actual uniform/graded PEC/dielectric operator controls agree with the basic
  engine within `2.8e-7` peak-normalized error after 120 steps. Two HIP domains
  and mixed CPU/HIP domains reproduce the new scalar CPU fields exactly in
  these controls. A deliberately missing halo produces error about `0.68`.
- Open MPI 4.1.6 tests pass with one/two/three ranks and multiple domains per
  rank. Two-rank GPU/GPU and separate CPU-only/GPU executables also pass over
  forced TCP loopback. A wrong epoch aborts the group with the expected error.
  Multiple physical hosts and AMD hardware have not been qualified.
- On two RTX 3080 devices, the static 25,165,824-node benchmark measured median
  5,752 Mcells/s on one device and 9,520 Mcells/s on two (1.66x), including
  staged exchanges. Resident field/coefficient arrays are 1,814,298,624 bytes
  for one domain and 1,816,657,920 bytes total for two equal domains. Driver,
  staging and operator-construction peaks are excluded from these byte counts.

These fixed-step controls establish the initial decomposition path only.
UPML, driven sources, opaque copper and production observations still need
their own cross-cut qualification before any full-antenna performance claim.
