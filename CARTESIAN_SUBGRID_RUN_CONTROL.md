# Native Cartesian-subgrid runtime control

Experimental native subgrids reuse the existing startup, runtime-pause and DFT
snapshot file protocol. This is an in-memory pause, not a restart checkpoint.
The owner drains pending processing and fences every logical worker at the
same completed E/H epoch before publishing `state="paused"`. Witnesses and
snapshot manifests include `nativeCompletedEpoch` and `nativeCompletedWorkers`.
The latter counts logical workers, not necessarily physical GPUs. Ordinary
engines report zero native workers.

Requests return after atomic submission. Controllers must wait for the matching
witness or `snapshot.xml`, and also check process liveness. A startup witness
can remain `paused` after a startup rejection or timeout has terminated the
process; the file alone is neither a liveness nor a convergence certificate.

```sh
openEMS input.xml --engine=gpu --exact-endcriteria \
  --runtime-control --decay-gate-id=example --decay-gate-timeout=600 \
  --dft-snapshots --decay-gates=-10,-20,-30,-40,-50
```

Use one controller and a fresh working directory on a filesystem supporting
atomic publication and hard links. The Python helpers are
`openEMS.run_control.control_run`, `openEMS.snapshots.request_dft_snapshot` and
`openEMS.startup_gate.decide_startup_gate`. A manual pause ID starts with the
launch's run token and a hyphen; continue/stop decisions must match the paused
witness ID. Existing Matlab startup/snapshot helpers are supplied; interpreter
execution remains separately unqualified.

Fields, material/PML histories, ownership, and DFT sample counts survive a pause.
Snapshots export accumulated sums without resetting them. Same-epoch requests
may therefore have identical values. Several decay thresholds can also be
crossed at one processing boundary. Runtime timeout resumes computation;
startup timeout/rejection fails qualification. A runtime stop ends gracefully
and is not energy convergence. Used snapshot IDs never overwrite published
data. Failed staging leaves diagnostics and no partial published snapshot;
a fresh request may succeed afterwards.

## Qualification and reproduction

The isolated Linux/CUDA R61 build was compared against the exact R60 runtime.
The full CPU suite passes76/76, all50 scripted controls pass, and16 frozen
operator/field arrays remain exact. Five native/ordinary protocol modes pass
on a3180-cell public fixture, retaining4000 steps, source, ports and timestep.
It tests control semantics, not radiation accuracy.

The unchanged152204-cell public dipole passes on one RTX4090 and an actual
RTX4090/RTX3080 pair. Each lane's160 raw final arrays are identical with and
without pauses and identical to accepted R59 arrays. Repeated snapshot sample
metadata is unchanged. Separate timeout/stop/write-error controls pass ten
child cases each on native CPU, two logical reference workers, and the actual
physical pair. Runtime-timeout and snapshot-error recovery preserve all final
arrays exactly. These are controlled identical-work comparisons, not claims of
bitwise equivalence between different backends or meshes.

The source and public fixtures are provided under `tests/subgrid_runtime/`.
With an already built matching executable, installed runtime/dependencies,
NumPy and h5py available in the test environment:

```sh
export SUBGRID_TEST_SOURCE=/absolute/path/to/openEMS
python3 tests/subgrid_runtime/application_control.py \
  /absolute/path/to/install/bin/openEMS \
  tests/subgrid_runtime/public-pause-fixture.xml /tmp/new-protocol-run reference-two
python3 tests/subgrid_runtime/application_fault_controls.py \
  /absolute/path/to/install/bin/openEMS \
  tests/subgrid_runtime/public-pause-fixture.xml /tmp/new-fault-run reference-two
```

For physical qualification, admit and select the two actual devices externally,
use `public-dipole.xml` and `cuda-two`. No launcher, device admission, remote
address, private antenna geometry, or credentials are embedded in these tests.
They start only their own solver children and terminate those children on a
test failure. The full public fixture has a larger setup/runtime cost than the
small protocol fixture.

Recent Metal, application-level backend-fault injection, opt-in overhead,
representative PCB memory scaling, curved-conductor spatial convergence and
full-array qualification remain open. R60's core wrong/missing/failed worker
acknowledgement controls do not replace application-level fault injection.
No production or frozen antenna runtime is changed by this increment.
