# LossyMetal setup diagnostics

Set `OPENEMS_PROFILE_LOSSY_SETUP=1` when running an openEMS build containing
PR242. The extension emits one `LossyMetal setup profile:` JSON record when
its setup returns. This is host setup timing; no device events, synchronization
or time-stepping instrumentation are added. The default is disabled.

```sh
OPENEMS_PROFILE_SETUP=1 OPENEMS_PROFILE_LOSSY_SETUP=1 \
  openEMS model.xml --engine=basic --numThreads=4 --no-simulation
```

`--no-simulation` prepares the operator only; its CLI exit code is 1 on the
successful preprocessing-only path. Check the completion record as well as
the exit code. No fields are advanced. In Python, the environment setting
applies to `FDTD.Run(path, setup_only=True)`; in Matlab/Octave it applies to
`RunOpenEMS(path, file, '--no-simulation')`.

Five disjoint phase durations cover surface scanning/classification, ADE
allocation, ADE coefficients, the interior-H mask and fit reporting. The first
phase includes initial state allocation and cached-line cleanup. Query counts
are exact. Query wall time samples the first call and every 1024th subsequent
call of each category, except fit solves, which are all timed. The output gives
sampled time and an estimated **inclusive** total; nested query estimates must
not be added together or treated as exact phase durations. Deterministic
sampling can alias regular geometry; inspect sample counts and compare coarse
phase times before drawing conclusions.

Disabled mode performs no profiler clock reads or logging. Enabled mode adds
counter increments and sampled clocks on setup paths; measure its overhead on
representative geometry. This instrumentation changes no conductor equations,
fits, mesh, materials or device update code. It does not by itself establish
that a performance problem is a numerical defect.

## Small public fixture

`benchmark/fixtures/lossymetal-cavity-1m.xml` is a preprocessing-only derivative
of the public TE101 cavity geometry in PR242's `python/Tests/Lossy_Metal.py`:
136×64×120 mesh lines, 1,012,095 geometric cells (1,044,480 grid locations), six
1-mm wall boxes. It uses smooth58MS/m copper, Rq0 and a Gaussian excitation
frequency definition. It has no source primitive and is not the original
cavity-Q validation. Run the command above in a disposable working directory,
using this XML path. No project antenna or supplier geometry is included.

On an Apple M5 Pro, four setup threads, three interleaved repeats against the
unchanged parent build: median setup0.4827s for parent,0.5153s for the new build
with profiling off,0.5059s with profiling on. The inactive instrumentation adds
about33ms (6.8%) on this short fixture; this is not a zero-overhead claim.
The on/off timing difference is within the observed short-run variability.
The change is diagnostic, not a setup acceleration or a production-runtime
recommendation. `--debug-operator --debug-PEC` outputs were byte-identical for
all three variants;9/9 CTest checks passed. No FDTD steps were run.

Recorded timings refer to their archived inputs and builds. Current public
fixture parameters are independent generic defaults; old timing records do
not establish new measurements for the revised fixtures.
