# Spatial Cartesian-subgrid curl blocks

Generated with Codex at the contributor's request. This is a partial C01/C03
increment in draft PR29, based on `d026d92`; it does not qualify a production
antenna, distributed material histories or larger cell guards.

## Representation and arithmetic

Ordinary four-edge face curls use spatial blocks, per-axis dimensions and a
mapping from structured E slots to the existing canonical independent edges.
Each leaf records its owning block and cell slot. A six-face mask identifies
regular rows; interface/PEC exceptions retain exact sparse patterns. Host
admission requires bit-identical ordered double coefficients and edge IDs.
No contour, coordinate, field order, interface mass or timestep is changed.
Borrowed metric views are supplied per call, so grid copies/moves do not retain
another grid's storage. Default CPU updates keep the previous fast pattern loop.
The native GPU preparation path selects the spatial host layout.

Device interiors resolve edge addresses from block coordinates and share exact
prepared FP32 high/low coefficient-only patterns. Curl and transpose read the
same pair. The transpose preserves canonical H order; exceptional rows preserve
the old sparse arithmetic. There is no coefficient rounding, independent
transpose fit, device double requirement, atomics or new backend. The portable
point operation serves reference/HIP/Metal. A per-case byte comparison retains
the previous device patterns when smaller. Rebuild C++ callers and backend
plugins against matching headers; mixed old/new plugins are not qualified.

An earlier unpublished candidate derived coefficients from binary-normalized
length/mass/reluctance factors. A broad scalar probe had small errors, but the
20,000-step field control failed the unchanged 2e-6 FP32-oracle gate on both the
reference backend and physical Metal. The rejected source/logs are preserved
locally. This increment uses exact prepared pairs instead; that failure is not
an established ordinary upstream solver defect.

## Focused evidence and reproduction

- `tests/test_cartesian_subgrid_regular_curl.cpp`: independent six-face matrix,
  one-ULP exception, failed-rebuild state preservation, exact expanded device
  operators, 128 driven CPU steps/copy/move, graded/material/multiple regions,
  mixed16/8/2 refinement and 64 regions.
- `tests/test_cartesian_subgrid_device_geometry.cpp`: malformed metadata and
  changed coefficient/transpose guards, seven 2,000-step fixtures and an actual
  structured graded-grid 20,000-step control. Fixed gates remain2e-4 against
  double CPU fields,2e-6 against the independent FP32 oracle and5e-4 undriven
  invariant drift. The structured path must actually be selected.
- Focused reference/physical Metal and separate physical RTX4090/RTX3080
  geometry controls pass. In the local 20,000-step case, fields match the
  independent FP32 oracle exactly and invariant drift is1.50e-8 or less.
- In that small graded case, curl/transpose metadata falls from620,896 to531,084
  bytes (14.47%); whole device allocation is597,124 bytes for the undriven case.
  Uniform fixtures retain the smaller old patterns. These are accounted bytes,
  not peak setup RSS, production capacity or a speedup claim.
- The earlier host-only r02 matrix passes13 stages, CTest38/38, Python14/14 and
  11,255 bit-identical prior arrays. A separate eight-case comparison retains
  all16 expanded-operator/CPU-field files exactly. This does not qualify later
  device changes automatically.
- The composed r04 local replay passes14 terminal stages, CTest39/39 and
  Python14/14;11,255 prior arrays are bit-identical. The full physical CUDA
  replay also passes33 gates and CTest39/39, with separate RTX4090/RTX3080
  geometry/material controls and retained lossless pair controls in both orders.
- Three interleaved same-build Metal repetitions compare equal duration and
  sampled final state. The spatial layout is slower: median1.28x at3,240 leaves
  and1.40x at25,920 leaves. It is a memory/addressing increment, not a throughput
  optimization. These are layout diagnostics, not an ordinary-master benchmark
  or antenna speedup. A separate unpublished transpose-slot candidate is being
  evaluated; it is not yet the implementation represented here.

With a configured test build:

```sh
cmake --build build -j4
ctest --test-dir build --output-on-failure -R 'regular_curl|device_geometry'
build/tests/test_cartesian_subgrid_device_geometry --metal
# A separately built HIP/NVIDIA lane, run once per actual device:
build/tests/test_cartesian_subgrid_device_geometry --cuda
```

Actual local r04 core SHA256:
`0acee95bbe7172a425207568c10ec6698c08a42fa9d682c35cf83c3ddf8d2db3`.
Actual r04 CUDA core SHA256:
`1a6fed18794cc748e148b10beeecbc88d14f771696f4247e9e8eafc928c1f7fa`.

Raw source/build manifests, negative attempts and numerical logs remain in the
contributor's evidence archive. Build hashes identify tested binaries, not a
claim that every caller/platform is qualified.

## Remaining work

Expanded CSR upload still exists as a bounded intermediate and independent
oracle. Host topology, fields, magnetic quadrature and observations retain their
reference organization. Device DFT, proportional multi-domain preparation,
distributed loss/PML histories, full structured scheduling, representative PCB
RF convergence, end-to-end memory/time savings and general interface scaling
remain open. No million-cell guard is raised. The production/frozen runtimes
and input models are unchanged.
