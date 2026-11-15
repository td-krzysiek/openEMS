# LossyMetal empty-space setup scan

This opt-in CPU preparation change is stacked on the worker/extent-cache
optimization. It does not change time stepping, surface-impedance fitting,
materials, grid coordinates or the ADE representation.

```sh
OPENEMS_LOSSY_SETUP=parallel OPENEMS_LOSSY_SKIP_EMPTY=1 \
  openEMS model.xml --numThreads=16 --no-simulation -v
```

`OPENEMS_LOSSY_SKIP_EMPTY` defaults off. With `1`, a Cartesian scan tests blocks
of up to eight Z lines before classifying individual E edges. The envelope
includes the staggered edge coordinates and is expanded by drawing tolerance
and outward floating-point rounding. A block is skipped only if every candidate
LossyMetal primitive is proven outside it. All three classification entries are
explicitly zeroed for each skipped node; subsequent magnetic-edge masks and
priority resolution retain their original behavior. A higher-priority void
can remove metal, but cannot create metal outside these candidate envelopes.
Zero-conductivity LossyMetal primitives are still candidates.

Unknown, transformed or unsupported bounds retain the full scan. Non-Cartesian
meshes and invalid drawing tolerances disable pruning. CSXCAD's optional
`CSPrimPolygon::GetConservativeBoundBox` API supplies cached Cartesian polygon
and linear-extrusion envelopes; these are prepared once, before workers run.
With older CSXCAD headers the code compiles and retains the conservative legacy
fallback, so polygon-heavy inputs may see little benefit. The legacy bounding-box
query contracts are unchanged. Rebuild dependent binaries when also adopting
CSXCAD's const-reference primitive-list API.

The setup log records `skip_empty`, `skipped_nodes` and `checked_blocks` once.
The state array and ordinary ADE storage remain allocated: this is not a sparse
field engine or a reduction of the simulation domain.

## Qualification

Use the fixed inspector from the worker/cache branch, then run:

```sh
python3 tests/check_lossymetal_setup.py \
  build/tests/inspect_lossymetal_setup /tmp/lossy-empty-check \
  /path/to/fixed-inspector-linked-to-reference
```

This compares complete ADE positions/coefficients/flags/counters/fits, base
VV/VI/II/IV and PEC dumps with pruning off/on, serial and parallel workers.
Controls include thin linear polygons, drawing-tolerance block boundaries,
transforms, curves, priority voids, Rq>0, missing thickness, zero conductivity
and unchanged PEC. An intentionally excluded-surface control must fail safely.

The diagnostic explicitly initializes PEC/PMC boundary-layer sizes and captures
ADE inside its live extension's `BuildExtension`. The earlier helper omitted
that initialization, then could access an inactive extension after the operator
had deleted it. Its empty-state receipts are not accepted qualification.
The diagnostic rejects PML; use the ordinary launcher for PML/FDTD controls.

Local qualification additionally uses two approximately one-million-location
fixtures, 4096-step basic/CPU/Metal controls with PEC and PML, and the standard
interleaved `benchmark/lossymetal_cpu.py` adapter. Full-antenna scaling and
CUDA/HIP execution remain separate checks. Frozen production builds are not
replaced by this change.
