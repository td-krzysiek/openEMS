# R60 — complete-epoch collective observation fence

This candidate adds an owner-thread fence before reusing the existing runtime
pause/snapshot protocol. It does not introduce another control-file protocol,
resume/checkpoint restart, or a field-migration/balancing mechanism.

`FenceObservationEpoch(expected)` rejects a half-step or an unexpected epoch.
Each local worker selects its actual backend and waits for completion of the
ordered device stream, including material/PML and DFT work. The coordinator
returns a logical-worker count only after all workers acknowledge that exact
epoch. This count is not a physical-GPU count. A missing, mixed or failed worker
acknowledgement poisons the coordinator; it cannot resume partial state.
The native operator exposes the same owner-thread boundary and invalidates its
observations/advance after an unexpected backend failure. Invalid caller epochs
are rejected before device work and do not poison a valid run.

No field, history or sum is downloaded by the fence. No barrier is added to the
ordinary stepping loop. This method blocks its owner thread; an external pause
request remains asynchronous and must be acknowledged only after the application
stops dispatch, drains asynchronous processing and successfully fences all
workers. That application-level integration is a subsequent increment. The CPU
native path is already synchronous. Stale H ghost copies at a completed epoch
are not changed: owned observations read owned state, and normal exchange occurs
before the next electric update.

## Predeclared qualification

- Full GPU-free CPU suite: expected 76 CTests, plus both focused fence targets
  and all retained parent controls, comparing 16 binary operator/field arrays
  against the actual frozen R58 prefix.
- Six direct-device cases: all three split axes, with/without source/resistors,
  one device, three logical workers, resumed 64-step fields and DFT. A forwarding
  test backend counts actual successful `Synchronize()` calls before the fence
  returns, avoiding a false pass caused by later blocking readback. Invalid
  half-step/epoch requests, wrong/missing worker acknowledgements and failures
  are negative controls. Existing nonzero domain error gates are retained;
  exact single-device equality is a same-runtime no-op control only.
- Twelve native cases: axes, shaped bulk loss/PML and one/two workers; 17 fences,
  64 resumed steps, exact matched E/H and complex DFT arrays. Counters must show
  no fence field download, DFT download or host-state allocation.
- Physical CUDA qualification is separate, on the previously authorized ordered
  and reversed RTX4090/RTX3080 pair and individual devices. It must follow exact
  accepted R58 physical and R59 radiation job identities. No overlap/admission
  inferred merely from a stale receipt.
- Metal remains separately qualified behind existing Mac headroom guards.

Source is a named R58 derivative. Production, antenna and frozen prefixes remain
unchanged. No claim of antenna convergence, external pause integration, speedup,
remote-worker transport or representative memory scaling follows from these tests.


## Scoped qualification — 17 November 2026

CPU76/76,16 exact actual-R58 operator/field arrays and80 physical invocations
pass.24 physical direct cases verify successful backend completion waits before
acknowledgement;48 native cases preserve resumed fields/DFT exactly, including
bulk losses/PML. Both orders of the actual GPU2RTX4090/GPU6RTX3080 pair pass.
The earlier GPU5 attempt held at the VRAM guard before docker run. Independent
review rehashes799 source/449 raw archived files; no solver source or gate was
changed for this resource-only pair derivative. External application pause and
snapshot integration is separate R61 and is not accepted from these core tests.
Recent Metal, scaling and representative antenna acceptance remain open.
Project evidence: `collective-fence-r60/independent-qualification-r01.json`.
