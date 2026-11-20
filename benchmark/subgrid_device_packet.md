# Explicit-device oblique packet replay

`subgrid_device_packet.cpp` is the byte-exact R72 device diagnostic used for
the accepted FP32-reference and physical Metal comparisons. It links against
a matching experimental openEMS source/build pair; no solver update is added.
Use the same dependency include/link flags as that build's control binaries.
For example, with the existing dependencies installed:

```sh
c++ -std=c++17 -O3 -ffp-contract=off \
  -I"$SUBGRID_SOLVER_SRC" -I"$SUBGRID_SOLVER_BUILD" \
  -I"$CSXCAD_PREFIX/include/CSXCAD" -I"$FPARSER_PREFIX/include" \
  benchmark/subgrid_device_packet.cpp \
  -L"$SUBGRID_SOLVER_BUILD" -lopenEMS -L"$CSXCAD_PREFIX/lib" -lCSXCAD \
  -Wl,-rpath,"$SUBGRID_SOLVER_BUILD" -o /tmp/subgrid-device-packet
PACKET_OUTPUT=$(mktemp -d)
/tmp/subgrid-device-packet --sampling-controls
/tmp/subgrid-device-packet 48 30 TE 2 2 2 .25 1.5 .7 "$PACKET_OUTPUT/reference.csv" --reference
/tmp/subgrid-device-packet 48 30 TE 2 2 2 .25 1.5 .7 "$PACKET_OUTPUT/metal.csv" --metal
```

Additional dependency include directories may be needed for a non-system
installation. Record the actual source, binary, loaded libraries and device.
The target is mandatory: `--reference`, `--metal` or `--cuda`; unavailable or
mismatched backends fail instead of falling back. CUDA is a prepared path,
**not qualified by the results below**. Admit memory and devices externally;
the standalone diagnostic is not a resource controller.

Arguments are base cells/axis, XY angle, TE/TM, independent refinement ratios,
longitudinal envelope width, wavelength, timestep fraction and CSV destination,
then the explicit target. The domain, seed, closed boundaries, receivers and
point/volume definitions follow [the CPU diagnostic](subgrid_oblique_packet.md).
This version retains its later 80 ps window and all 15 history columns, with
separate electric/magnetic timestamps. It reads back every complete epoch;
there is no throughput or production-observation overhead claim.

Eight accepted cases use 48 base cells per axis and wavelength1.5:

| Polarization | XY angle | X/Y/Z ratios | Envelope width | Timestep fraction |
| --- | ---: | --- | ---: | ---: |
| TE | 30 degrees | 2/2/2 | .25 | .7 |
| TM | 30 degrees | 2/2/2 | .25 | .7 |
| TE | 45 degrees | 2/2/2 | .25 | .7 |
| TM | 45 degrees | 2/4/1 | .25 | .7 |
| TE | 30 degrees | 2/2/2 | .25 | .35 |
| TM | 30 degrees | 2/2/2 | .25 | .35 |
| TE | 30 degrees | 2/2/2 | .18 | .7 |
| TM | 45 degrees | 2/4/1 | .18 | .7 |

On one Apple M5 Pro, independent review of every raw point/volume history gives
maximum physical-Metal differences of5.152e-7 from the frozen double CPU and
4.484e-7 from the FP32 reference, below the respective2e-4 and2e-6 gates.
All original resolved volume5%, temporal1% and invariant5e-4 gates pass.
Maximum child RSS is597,032,960bytes. A reviewer initially selected the CPU
owner file instead of the pinned Metal owner; that failure is preserved and
the corrected review retains the exact original numerical/resource comparison
block. No solver or raw histories were changed or rerun for that correction.

These are backend-agreement controls. Original point-H spatial accuracy
failures remain failed; agreement between backends cannot repair them. The
finite closed-box tail does not separate interface reflection or establish
pulse completion. Other incidence planes, an independent ordinary engine,
broader interface/feature qualification, physical multi-device and representative
antenna validation remain open.
