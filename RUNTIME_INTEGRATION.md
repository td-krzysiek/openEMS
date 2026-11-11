# Upstream / multi-GPU integration — 11 November 2026

Generated with OpenAI Codex at the user's request. This is an experimental fork
integration, not a deployment or an upstream submission.

The branch combines upstream `1d710c210920908a143b573f2babca05455f5479` with
our adopted runtime `836dd5d29b7afa3f2b96b3a5b2ca6eef2712b832`. It retains
the upstream operator/coefficient refactor and the experimental CPU engine.

## Retained extensions

- Native HIP and Metal X-domain execution, fixed uneven cuts, owned energy
  reductions and post-extension E/H halo exchanges.
- Per-domain Drude/Lorentz and PR242 LossyMetal material histories.
- Experimental opaque conductor boundary, shared coefficients and compact GPU
  histories, plus compact HIP UPML coefficients.
- Startup qualification, manual/decay pause gates and accumulated DFT snapshots.
- Optional HIP/CUDA NF2FF and Python/Matlab controls.

## Integration choices

The old experimental `--engine=tiling` and its extension hooks are removed.
The upstream `--engine=cpu` provides the maintained temporal-blocking path.
The opaque extension runs flat on that CPU engine; it does not claim support
for its box hooks. Its HIP/Metal device updates remain available.

The GPU operator no longer has a dedicated subclass or full intermediate
coefficient arrays. Domain uploads read the upstream coefficient interface.
Domain exchanges still finish after the device extension hooks. Upstream
single-GPU cylindrical sub-grid stepping is preserved; the domain adapter
continues to reject cylindrical meshes and unsupported device extensions.

Rebuild Python bindings and GPU plugins together with the solver. Use a new
prefix: do not replace a runtime referenced by frozen or active simulations.

## Qualification

Build, numerical and runtime-control checks are required before adoption:

1. CTest and Python API/unit regressions.
2. Single versus two/three domains with PEC, PML, RLC, opaque and dispersive
   materials, including PR242. Check all nonzero complex E/H components.
3. Uneven cuts, one-node slabs, invalid partitions and unsupported extensions.
4. Startup/runtime pause, DFT snapshots, continuation and energy stopping.
5. Old/new public antenna input and native CUDA tests on free physical GPUs.

The Metal test plugin exposes logical domains on **one physical GPU**. It is
not installed and cannot establish physical multi-GPU scaling. Preserve failed
controls and distinguish inherited upstream behaviour from integration defects.
Only measured setup/FDTD times may be reported as performance evidence.
