# Experimental engine regression

Use a separately built original and candidate library with matching Python
bindings and dependencies. Do not install over your production solver.
Run from the repository root, with the installed Python environment active.

```sh
LD_LIBRARY_PATH="$BASELINE_LIB:$PREFIX/lib" python TESTSUITE/performance/validate_engine.py run --case dispersive --threads 4 --output /tmp/openems-reference
LD_LIBRARY_PATH="$CANDIDATE_LIB:$PREFIX/lib" python TESTSUITE/performance/validate_engine.py run --case dispersive --threads 16 --output /tmp/openems-candidate
python TESTSUITE/performance/validate_engine.py compare --reference /tmp/openems-reference --output /tmp/openems-candidate
```

Repeat for `pml` and `mixed`, with 4 and 16 threads and fresh output paths.
Tolerance per array: maximum absolute error <= 1e-6 * reference peak + 1e-20.
The fixture rejects non-finite arrays and trivial excitation. These short
2500-step regressions are not antenna convergence or performance benchmarks.

Run the upstream unit suite and Python smoke test against the selected library:

```sh
python -m unittest discover -s python/Tests -p 'test_*.py' -v
python .github/smoketests/python/MSL_NotchFilter.py
```
