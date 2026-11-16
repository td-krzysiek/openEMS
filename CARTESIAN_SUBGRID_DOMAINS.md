# Experimental compact-domain execution

This is a scoped C++ core interface. The
[canonical checklist](CARTESIAN_SUBGRID_CHECKLIST.md) separates implementation
from physical backend, processing and representative-antenna qualification.
The native XML/Python/Matlab engine selector remains single-device.

## Ownership, workers and transport

`CartesianSubgridDomainPlan` assigns fixed ownership from one axis and a list
of increasing internal coordinate cuts. Cuts assign work; they do not create
mesh faces. Up to 64 nonempty domains are admitted. All six magnetic
quadratures of a leaf stay together. Each coupled electric mass component
belongs entirely to one worker, so geometric assignment may move electric
edges across a cut. The plan records that count. This avoids splitting a
coupled solve; it does not guarantee balanced work.

Each worker has owned-first local E/H arrays followed by read-only ghosts.
Local row remapping retains the exact operand order and split coefficients.
Ghost dependencies may fan out or cross non-neighbour workers. Only unique
owned export values are gathered, downloaded and copied into preallocated
receive buffers, then contiguous ghost ranges are uploaded. No whole-field
transfer is performed at every half-step. The plan still uses a temporary
global expanded upload during setup; production-proportional setup memory
and structured interiors remain open.

The coordinator knows only worker and transport interfaces. GPU handles and
local gathers/uploads belong in the local worker adapter, not the scheduling
logic. A reference worker runs on the CPU; this is not integration with the
optimized production CPU engine. Remote workers, MPI/RPC and asynchronous
transports are future adapters, not implemented capabilities.

At epoch n, exchange H before the E update; then exchange E at n+1 before the
H update. Workers reject missing/wrong-phase ghost acknowledgements. Full-step
observations read owned fields and scatter each global state exactly once.
Queued epochs are logical submission state, not device fences; export/owned
reads and ghost uploads finish before returning. Errors after a submission
poison the coordinator so it cannot silently continue a partially changed step.
There is no single borrowed device-field buffer for the whole domain run.

Use `FDTD/cartesian_subgrid_domains.h` from matching source headers. For example,
`CreateCartesianSubgridDomains(plan, backends)` builds local workers and a
host-staged transport. The more explicit worker/transport/coordinator factories
permit alternate adapters without putting communication choices in the core.
Backends must outlive workers. The optional `DeviceUpload` backend operation
fails closed before worker allocation when unsupported; the existing one-device
factory remains usable. This changes C++ layouts/vtables: rebuild dependent
callers/plugins with matching headers. Built-in backend tests do not qualify
plugin mode or Windows symbol export.

## Qualification scope

The [domain regression](tests/test_cartesian_subgrid_device_domains.cpp) covers
one through four logical domains, XYZ cuts, separate regions, anisotropic
ratios, nonuniform meshes, a lossless dielectric step, PEC and soft sources.
It compares complete fields with the existing one-device FP32 execution and
double CPU reference over 1500 steps. It also checks actual transfer payloads,
stale/missing exchanges, poisoned continuation, invalid uploads and unchanged
one-device operation when upload capability is absent. Physical backend
requests reject fallback rather than accepting a reference result.

The immutable local suite passes all 19 stages: CTest 23/23, Python unit modules
14/14, the prior native CPU/reference/Metal fields, ports/DC, PML, dipole,
polygons and ordinary-engine controls. Thirteen previous core/NF2FF prefixes
remain unchanged; 11,255 numerical arrays equal the frozen larger-interface
build exactly. Seven reference cases and eight physical Metal cases pass,
including mixed CPU-reference/Metal workers. Across these 1500-step domain
controls, one-device fields are identical and maximum double-reference error
is 1.62522e-6 (gate 2e-4). Missing/stale exchanges and poisoned continuation
reject. A separate ASan/UBSan domain-source control passes; LeakSanitizer is
unsupported on this macOS and no leak-sanitizer claim is made.

A separate 360,424-leaf control exercises an actual 16,384-edge mass component
across two logical Metal domains on one M5 Pro: base 8194/4/4, spacing
0.001/0.0013/0.0009m, region (2,1,1)–(8192,3,3), ratio 2/2/2 and an X cut
at 4.097 m. After 256 steps at 0.8 of the stable timestep, the one/two-domain fields
are identical, double error is 1.02301e-6 and invariant drift 1.18914e-9.
Some 32,772 electric states move owner to preserve whole coupled components.
These are functional controls, not a scaling benchmark.

Frozen local prefix: `output/openems-subgrid-domains-qualified-20261115/install`.
Configure-time source parent is `1a8f48f`, with this domain increment. Exact
source inventories, logs and identities remain in the project evidence tree:
`experiments/openems-cartesian-local-subgrid-20261112/domains/`.

- Core SHA256: `ebb022b5539b336a6b4b43b06bc1e83744b3a657ebd77d8bbb3e10489648ba53`.
- NF2FF SHA256: `87890b7ed62795695b7f81e15e3453f3c789c4f75821e5508bb1157c869c27ef`.

The same numerical source passes eight logical-domain controls on one physical
RTX3080, including mixed CPU-reference/CUDA workers, through existing HIP/nvcc.
One-device errors are zero; maximum double error is 1.62522e-6. CTest 23/23 and
single-device fields/ports/DC pass. The first attempt hit an external 300 s
CPU-test deadline; the second retains identical source/runtime/gates, runs CUDA
first and allows 900 s for CTest. Both attempts remain recorded. The final container
exited 0 without OOM, and final source/runtime hashes match. This does not
qualify physical multi-GPU, domain ports/PML or native CAD CUDA paths.

- CUDA core SHA256: `de9994d3e5a09f07de74de1effed5cc9c6371d2d7ee28bd0204d8ccdeab56c0a`.
- [CUDA reproduction scope](CARTESIAN_SUBGRID_CUDA.md#scoped-logical-domains-on-one-rtx3080).

All existing total-cell, interface-factor and refinement-alignment guards
remain. The domain constructor rejects resistor ports/history and the local
upload rejects nonzero unsupported ledger state. PML, losses/conductors,
distributed ports, DFT/NF2FF, collective pause/resume and native processing
are not implemented in this domain lane. Their one-device tests do not imply
distributed support. No throughput, representative PCB accuracy or physical
multi-GPU claim follows from logical domains on one physical device.
