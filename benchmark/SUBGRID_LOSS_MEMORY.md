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

On one Apple M5 Pro, two AB/BA samples reduced median peak process RSS from
84.77 to 51.77 MiB for bulk and from 251.41 to 193.54 MiB for surface traces.
Printed samples/energy were identical and all 83 CPU/reference CTest cases
passed. These short measurements do not establish a runtime speedup, whole-field
bitwise equivalence, physical-device equivalence or full-antenna capacity.
The original full application remains held after a 64 GiB preparation OOM;
its other dense arrays and device-packing temporaries are separate work.
