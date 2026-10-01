# Experimental Cartesian tiling engine

This branch ports the experimental engine from
[biergaizi/openEMS, project-diamond-rework1](https://github.com/biergaizi/openEMS/tree/project-diamond-rework1)
(commit `289fef9`) onto upstream `7ae67b5`.
The original tiling implementation is by Yifeng Li; existing file copyright and
GPL notices are retained. This port and its tests were substantially AI-generated
with OpenAI Codex and require human review before upstream submission.

The standard engine remains the default. Opt in with:

```python
fdtd.Run(sim_path, engine='tiling', numThreads=4)
```

```matlab
RunOpenEMS(sim_path, sim_file, '--engine=tiling --numThreads=4');
```

```bash
OPENEMS_TILE_X=10 OPENEMS_TILE_Y=10 OPENEMS_TILE_STEPS=5 \
  /path/to/candidate/bin/openEMS model.xml --engine=tiling --numThreads=4
```

`OPENEMS_TILE_X` and `OPENEMS_TILE_Y` are minimum diamond widths in grid cells
(default 10, range 1–4096). `OPENEMS_TILE_STEPS` is the number of full E/H steps
per temporal block (default 5, range 1–64). The actual diamond footprint grows
across its time window. Z is never partitioned because the SIMD layout packs
that direction. Remainders of each requested iteration batch use synchronized
single-step spatial tiling; probes/dumps are sampled at their original times.

## Supported scope

Cartesian grids; ordinary PEC/PMC boundaries; UPML; soft electric/magnetic
excitation; lossy nondispersive materials; Lorentz/Drude and conducting sheets.
The local numerical fixtures cover UPML, electric port excitation, lossy
substrate, conducting sheet, and two-order electric/magnetic dispersion.
Do not infer coverage of every excitation, boundary or material combination.

Unsupported active extensions (including Mur, TFSF, reactive lumped RLC,
steady-state detection) fail during setup before workers start. Ordinary lumped
resistive ports work; an empty RLC extension is allowed. Cylindrical coordinates
are rejected. No MPI integration is provided. Do not use this branch for
production antenna results without a matched reference/convergence check.

## Review changes relative to the historical branch

- Preserve current compressed-array storage and current material equations.
- Keep tiling opt-in; retain existing engine enum values.
- Keep Z whole and distribute work over XY, avoiding the historical choice of
  an unpartitioned Z axis for worker allocation.
- Retain E-only or H-only clipped edge updates. The old combined planner could
  discard both when only one half of a tile was empty.
- Initialize extension maps before worker startup; reference index lists rather
  than copying vectors in each update. Skip duplicate range maps.
- Validate tile settings, report setup errors, and reject unsupported physics.
- Shut down workers through their start barrier. A startup gate also releases
  partially created teams if thread construction fails.
- Remove unused demonstration planners/debug output and nonstandard variable
  length arrays from the imported planner; retain C++11 and existing dependencies.

The baseline update functions remain separate from the tiled range functions.
This intentionally avoids changing the default engine's hot loop during the
experimental port. Their equations must be kept aligned in future changes.

## Local tests

Run outside `python/`, using the matching installed Python interface and candidate
library path. Keep reference and candidate installations separate.

```bash
c++ -std=c++11 -Wall -Wextra -pedantic -I. \
  TESTSUITE/tiling/test_plan.cpp tools/tiling.cpp -o /tmp/tiling-plan
/tmp/tiling-plan

export LD_LIBRARY_PATH=/path/to/candidate/lib:/path/to/dependencies/lib
python3 -m unittest discover -s python/Tests -p 'test_*.py' -v

# Optional retained numerical evidence (use a new directory each time).
OPENEMS_TILING_TEST_OUTPUT=/path/to/new-results \
  python3 python/Tests/test_tiling.py -v
```

The planner test checks exact E/H coverage for 2,280 combinations, including
one-cell tiles, partial edge tiles and blocks larger than the domain. It can be
compiled with AddressSanitizer/UndefinedBehaviorSanitizer.

The Python regression launches separate bounded processes and compares port
traces, port spectra and E/H field dumps. Acceptance is peak-normalized error
`<= 1e-6 * reference_peak + 1e-20`, with shape and finite-value checks. It tests
1/2/4 workers, one-step/temporal blocks, remainder batches and explicit Reset.
It also checks clean rejection/destruction for unsupported physics and invalid
settings. These tests are numerical regressions, not RF mesh convergence tests.

## Remaining acceptance gates

1. Repeat regressions on the target native ARM Linux compiler/runtime.
2. Run sequential, alternating A/B fixed-step 2x2 benchmarks with identical XML,
   CPU affinity, compiler flags, dumps and memory policy; measure setup time and
   peak RAM as well as FDTD throughput. Test 8/16 workers and 2/5/8-step blocks.
3. Compare a complete converged antenna run: port traces/S11, convergence history,
   far-field data, stopping schedule and wall time. Never label truncated runs
   as antenna results.
4. Run Octave integration tests and thread-race checks; x86/SSE validation remains
   required before claiming portable equivalence. The partial thread-creation
   failure path still needs fault-injection coverage.

No production speedup is claimed by this port.
