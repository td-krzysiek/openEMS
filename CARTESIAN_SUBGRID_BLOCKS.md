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
  or antenna speedup. The subsequent exact-order and streamed-preparation
  increments are recorded below.

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

## Exact ordering and streamed preparation

The r06 increment stores an exact0..3 edge-slot tag in regular transpose
terms and the canonical four-term order in the existing face descriptor.
It removes repeated edge searches and per-step sorting, retaining exact
coefficient pairs and sparse exceptions. Physical RTX4090/RTX3080 replay
passes33 gates/CTest39/39, including both retained lossless device orders.

The composed r07 single-device preparation streams canonical rows directly
into translated patterns. Transpose setup retains incoming face IDs, then
recovers the same coefficient pairs from the compressed curl. It no longer
retains expanded column/high/low arrays for both operators while packing
alternatives. The native conductivity/surface provider uses the state-only
path too. Expanded `DeviceData()` remains an independent arithmetic oracle
and is retained by domain preparation; unsupported mixed state/CSR requests
reject explicitly. Dictionaries are temporary; the timestep kernels and
byte-selection rule are unchanged from r06.

The composed r07 local replay passes14 stages, CTest39/39, Python14/14 and
11,255 unchanged prior numerical arrays. Tests compare selected packed
operators, exact coefficient pools, order and state-only factor/source/port
metadata against independent expanded preparation; malformed/incomplete row
streams and duplicate ordering are rejected. Physical Metal includes the
20,000-step exact-FP32 field and strong-DC controls. The separate r07 physical CUDA replay also passes33 gates and CTest39/39,
including each RTX4090/RTX3080 and both retained lossless pair orders. Actual
Docker terminal state is exit0/no OOM; source/runtime identities are rechecked.

Three interleaved Metal repetitions on a491,520-leaf procedural fixture give:

|Median|r06 expanded preparation|r07 streamed preparation|
|---|---:|---:|
|Peak process RSS|1,614,561,280B|1,233,518,592B|
|Whole setup|3.0630s|3.1654s|
|2,000 timesteps|2.3479s|2.3511s|
|Device allocation|141,076,208B|141,076,208B|

Peak RSS is23.60% lower; setup is3.34% longer. This is a memory improvement,
not a timestep speedup. The sampled final fields/energy, device allocations
and timestep are identical. A larger interleaved attempt was stopped by its
owned-process memory-reserve guard and remains partial evidence.
The r06 spatial kernel still takes about9.78% longer than old translated
patterns at699,840 leaves while using26.25% less accounted device storage.
No production/ordinary-Yee performance claim follows from these diagnostics.

Actual r07 local core SHA256:
`4af61e6916978ba5d017b802e1bc64355cd6b3ce133d3753568e4dde896ece1f`.
Actual r07 CUDA core SHA256:
`a69a976f05a9309a2746af23adbab5d0962a3bc4b657bff045ba82ae11d88d94`.
The native replay is a correctness suite, not a CPU speed benchmark; its
reference-backend timings need separate controlled review before performance
adoption. No allocation guard has been raised.

## Remaining work

Expanded CSR remains the independent oracle; default single-device, global
domain and local domain preparation now use streamed paths. See the
[local domain setup measurements](CARTESIAN_SUBGRID_DOMAINS.md#streamed-local-domain-preparation).
Host topology, fields, magnetic quadrature and observations retain their
reference organization. Device DFT, proportional multi-domain preparation,
production-scaled loss/PML histories, full structured scheduling, representative PCB
RF convergence, end-to-end memory/time savings and general interface scaling
remain open. No million-cell guard is raised. The production/frozen runtimes
and input models are unchanged.
