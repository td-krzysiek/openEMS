# Experimental subgrid loss-preparation memory

Build with `OPENEMS_BUILD_TESTS=ON`, then:

```sh
cmake --build build --target benchmark_subgrid_loss_memory
OPENEMS_SUBGRID_SETUP_PROFILE=1 /usr/bin/time -l build/tests/benchmark_subgrid_loss_memory bulk 40
OPENEMS_SUBGRID_SETUP_PROFILE=1 /usr/bin/time -l build/tests/benchmark_subgrid_loss_memory surface 250000
```

Use `/usr/bin/time -v` on Linux. Run parent/candidate in both orders on the
same machine and compare maximum resident set size, not just allocator payload.
The bulk case has 64,007 leaves and constant electric conductivity. The surface
case has 71 leaves and repeats 250,000 passive seven-pole paths; it stresses
representation capacity, not a physical metal fit. Both prepare the actual
operator and take sixteen CPU steps. Their printed energy/field samples are
diagnostics; the CTest surface/volume oracles provide independent numerical
coverage. No antenna input or geometry is required.

For R83 on one Apple M5 Pro, two AB/BA samples reduced median peak process RSS from
84.77 to 51.77 MiB for bulk and from 251.41 to 193.54 MiB for surface traces.
Printed samples/energy were identical and all 83 CPU/reference CTest cases
passed. These short measurements do not establish a runtime speedup, whole-field
bitwise equivalence, physical-device equivalence or full-antenna capacity.
The original full application remains held after a 64 GiB preparation OOM;
its other dense arrays and device-packing temporaries are separate work.

## Device loss-table assembly

The optional `--pack` mode prepares cold device/domain loss tables and writes
`loss-tables.bin` in the working directory. The diagnostic encoding contains
the timestep followed by native length-prefixed U32/F32/F64 arrays; it is not
a portable simulation checkpoint. Use a fresh directory for each run:

```sh
bench="$(pwd)/build/tests/benchmark_subgrid_loss_memory"
case_dir="$(mktemp -d)"
(cd "$case_dir" && "$bench" sparse 48 --pack)
```

`sparse 48` uses 110,599 leaves and 25,000 passive paths, leaving most electric
rows empty. Its repeated generic voltage-path resolution dominates the whole
command time; this is not a native-metal setup speed benchmark.

R84 compares the R83 parent with counted CSR assembly, using two AB/BA samples
per case on the same Mac. All emitted loss arrays agree byte-for-byte; physical
field/energy diagnostics agree and all 84 CPU/reference CTest cases pass.

| Mode with `--pack` | Parent median peak RSS | Counted median peak RSS |
| --- | ---: | ---: |
| bulk48 | 149.41 MiB | 115.55 MiB |
| surface250000 | 404.85 MiB | 294.02 MiB |
| sparse48 | 111.52 MiB | 93.93 MiB |

These are process RSS measurements, not device VRAM or full-application
capacity. Actual Metal, CUDA/multi-device and representative antenna gates
remain separate. No timestep kernel was changed.
