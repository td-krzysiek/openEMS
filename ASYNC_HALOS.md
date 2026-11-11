# Experimental asynchronous host-staged halos

New, separately qualified HIP domain runs can select:

```sh
OPENEMS_DOMAIN_HALOS=overlap openEMS model.xml --engine=gpu --gpu-devices=0,1
```

The default is `serial`. `async` uses reusable pinned staging, copy streams and
completion events while retaining full half-step barriers. `overlap` additionally
submits the H interior while the corrected E halo is in flight, then completes
the H boundary after arrival. Each H node is updated once, with the same curl
arithmetic and timestep. It is an E-to-H pipeline only; there is no speculative
next timestep or dynamic repartition.

Both field phases still exchange tangential neighbor planes every timestep.
The full-step barrier remains before processing, DFT, energy, pause and snapshot
handling. Extension order is unchanged. The experimental native domain backend
retains its existing allowlist of local extension implementations. Current
pre-update hooks must not consume an incoming E halo before `BeginPhase(Current)`.

The coordinator/executor interface contains no HIP API. The CPU domain fixture
implements the same disjoint current ranges for numerical/scheduling regression.
In overlap mode its Voltage observation happens after E arrival, while H interior
work may already be in flight. Observers requiring a consistent E/H time pair
must use the full-step boundary, as the native Engine processing already does.
Other device adapters explicitly reject unsupported opt-in modes. Existing
single-GPU and default Metal/reference paths retain their update scheme.

HIP allocates two-component pinned buffers, copy streams and reusable events at
setup. A producer event orders downloads after the native E/H hooks. All neighbor
downloads are submitted before any completion is waited for; all ready arrivals
are given a chance to enqueue their uploads before destination waits. Local
transport retains its host-buffer wire contract, including the host copies.
No P2P, peer access or extra full-field allocation is required. Five domains with
Y/Z counts1611/106 require about20.85MiB of additional pinned staging RAM.
Allocation/transfer failures are errors, not silent fallbacks.

Rebuild core and plugins together: the optional backend interface has changed.
A frozen running runtime must not be replaced. Profile headers record
`halo_mode`; exchange wall scopes include deferred arrival waiting, overlap
compute and are not pure DMA time. Pack/unpack scopes in async mode describe
submission/host copying, not synchronous transfer duration. Compare total FDTD
wall time, fields and probes before adopting a mode.

## Reproduce

Build/test the dependency-free coordinator:

```sh
cmake -S FDTD/distributed -B build-runtime -DOPENEMS_DOMAIN_TESTS=ON
cmake --build build-runtime
ctest --test-dir build-runtime --output-on-failure
```

With rebuilt HIP binaries and two exclusively reserved physical GPUs:

```sh
python tests/manual/async_halos.py --binary /prefix/bin/openEMS \
  --devices 0,1 --output /new/qualification --benchmark
OPENEMS_DOMAIN_HALOS=overlap python tests/manual/runtime_gates.py \
  --binary /prefix/bin/openEMS --engine gpu \
  --solver-option=--gpu-devices=0,1 --output /new/gates
```

Public fixtures cover vacuum/PML, PEC, RLC, opaque metal, Lorentz material,
LossyMetal and a patch. Compare complex fields/probes at1e-5 peak-relative and
patch reflection at1e-4 absolute; record actual errors. Thin ownership tests
include a one-line domain. Benchmarks rotate serial/async/overlap order on one
machine. Logical-device runs do not establish physical multi-GPU performance.
No speedup is claimed before physical-device qualification and measurement.
