# Antenna runtime integration

Snapshot recorded on 25 October 2026. This branch preserves the adopted
experimental solver changes and the qualified startup gate. It is not a merge
of every prototype branch and is not the latest upstream solver.

## Source identities

| Fork branch | Source | Status |
| --- | --- | --- |
| `runtime/monster-20261025` | `8a6afe4937bc4493a8bb32bb8c542411da909393` | Exact solver source of the active antenna workers when checked on 25 October |
| `integration/antenna-runtime` | `20635e910cfde2399b08d5085db600fbb86526c1`, followed by this documentation | Same integration plus opt-in startup qualification; qualified on short controls, prepared for future antenna runs |
| `upstream/sync-20261025` | `d686f1dfbb8900a748eacf1925e091ee32d1995c` | Unmodified upstream master snapshot; not merged or deployed |

The integration retains opaque-surface timestep stabilization, connected box
unions and CPU/Metal/HIP updates; shared material coefficients and compact GPU
history storage; compact HIP UPML coefficients; setup-error cleanup; accelerated
HIP NF2FF with the CUDA build path; and on-demand accumulated DFT snapshots.
Earlier experimental CPU tiling also remains in its ancestry. The optional
CUDA FMA experiment and multi-domain/MPI proof of concept are separate branches,
not additions to this runtime integration.

The startup gate is the only solver commit after the active runtime snapshot.
It continues the same process after qualification without rebuilding the
operator, resetting fields or discarding DFT sums. See [STARTUP_GATE.md](STARTUP_GATE.md)
and [fork PR 18](https://github.com/td-krzysiek/openEMS/pull/18).

## Build provenance and qualification

Source checkout alone is not a reproducible installation: CSXCAD, compiler,
build options, GPU architecture and linked libraries must also be pinned.
The observed active Linux runtime has these SHA-256 hashes:

```text
openEMS       b22dc28dd7caec8299de56fcf7f8c6f2cbc36e5f037e727413e0975d91fc062e
libopenEMS.so 15fd79d5f3589ae097710784583b2cd381d71fb604f10b2ffa81508e7eb4afac
```

The separately qualified startup-gate Linux build has:

```text
openEMS       02bc7d8c0e073bde65ad839b82983d015fc59e09dba1d4953c1ce89050f36f63
libopenEMS.so f1c49ddd1a3004b333d4c077b8c7cfbde43a8dcb31845bf33d19d0bb6d25f0ce
```

Startup continuation controls passed on basic CPU, GPU reference, actual Metal
and RTX 3080 CUDA, including opaque histories, rejection, timeout and abort.
Comparisons against uninterrupted controls used stated tolerances and measured
zero field/probe differences. A production-controller synthetic fixture also
matched the deployed runtime on 226 probes and 12 nonzero complex field arrays.
These are short regression controls, not a new full-array antenna qualification.
Matlab helpers were reviewed but not executed. Rebuild bindings and C++ clients
against the new library; do not mix ABI versions.

This backup does not change active workers or replace frozen antenna results.
It contains solver source and synthetic tests, not antenna designs, supplier
documents, credentials, compiled libraries or simulation result archives.

## Upstream review

The common ancestor with upstream is `c3a956b221ee57f85bcf45f65f853c06e9769a7d`.
The fetched upstream snapshot contains 19 commits beyond that ancestor,
including 12 beyond the earlier `11aa063` review baseline. Dates below follow
commit metadata; the comparison uses ancestry, not date ordering.

| Change | Commits | Relevance to this runtime |
| --- | --- | --- |
| HIP runtime plugins, unified GPU selection, early device checks | `170c71b`, `c0ca95a` | Cleaner deployment and earlier failure when no usable device exists. Changes backend ownership and construction; port opaque extensions and requalify selection. It does not split one solve across GPUs. |
| Builds without a local GPU, Windows plugins and GPU build CI | `81ec4db`, `a5594df`, `4ed0ced` | Useful for portable cloud/container images. Architecture compatibility and plugin build identities still need checking. |
| Compact HIP UPML storage | `53a4264` | Based on our upstream PR 238. Same purpose as our existing memory optimization; reconcile implementations rather than stack another copy. |
| GPU synchronization before returning from a timestep batch | `6bb955d` | Fixes inflated solver speed reports when queued work was not timed. Does not itself promise faster solves; remeasure throughput and wall time. |
| Binary benchmark harness | `6fa7e8e` | Five workload families at approximately 1M/8M cells, interleaved binary comparisons and probe checks. Useful for future performance PRs; does not replace antenna/opaque-boundary controls. |
| Mesh-derived geometry tolerance | `db68e62` | Tolerance becomes 1/1000 of the smallest mesh step. Recheck thin copper, contacts and mesh-aligned unions; geometry classification can change. |
| C++17 and VTK dependency removal | `a2b049d`, `d9ca3b6`, `11aa063` | Update build images and bindings. Internal VTK output uses zlib; this reduces dependencies, not necessarily FDTD runtime. |
| Safer setup, CLI and timestep limits | `dd72017`, `c971899`, `ec4d2db`, `7badca4` | Better failure handling and an INT_MAX timestep cap. Check worker option parsing and error paths. |
| Delayed excitation correction | `d1fef7c` | Skips excitation before its delay expires on CPU/HIP/Metal. Particularly relevant to delayed step/custom sources; inspect actual source settings before attributing antenna changes to it. |
| Optional DC-free Gaussian excitation | `d686f1d` | Addresses static residual charge for pulses with significant DC content. Opt-in; do not change frozen antenna excitation settings implicitly. |
| Setup-cleanup test timing | `ff68e3d` | Allows thread counts to settle before checking them; test robustness rather than a numerical change. |

Primary evidence is the corresponding commit diff in
[upstream openEMS](https://github.com/thliebig/openEMS/compare/c3a956b...d686f1d).
This review reads the source; it does not claim the new upstream was built or
numerically qualified against our extended runtime.

Before adopting a merged runtime, reconcile plugin ownership and UPML changes,
build with matching CSXCAD and bindings, then run CPU/Metal/CUDA controls for
opaque stability/unions, allocation cleanup, DFT snapshots and startup
continuation. Compare both numerical output and external wall-time benchmarks.
