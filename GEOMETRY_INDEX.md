# Experimental setup geometry index

## Prepare once, borrow throughout setup

With CSXCAD prepared-grid support, select `OPENEMS_GEOMETRY_INDEX=prepared`.
One operator-owned lookup prepares shared exact candidate lists before material
and PEC construction. Materials, material|metal and excitation retain separate
priority orderings inside that object. Each mask covers the three full-axis
query orientations, including unrestricted-axis sentinels. This covers built-in
line/plane queries used by material/PEC, UPML, dispersive/conducting-sheet
extensions, excitation, Mur, cylindrical setup and setup material dumps.

Built-in consumers call `GetPrimitivesBoundBoxView(scratch, x, y, z, type)`:
registered queries return immutable vectors by reference, without geometry
queries, merging, filtering, allocation or copying. Identical lists share
storage within each prepared grid. Diagnostics print total mapping bytes,
per-grid unique-list counts and primitive references. Unsupported masks,
fully bounded XYZ boxes, calls after setup, and non-prepared modes retain the
legacy virtual query through caller-owned scratch. The original by-value API
remains available; it does not provide zero-copy reads.

Custom operators overriding the original query must also override the new view
method if they enable prepared mode and need custom registered-query behaviour.
New virtual methods require rebuilding dependent native code. The table has the
same fixed-geometry/lifetime requirement as the spatial snapshot and is freed
before FDTD, including exceptional setup exits. This prepares candidate lists,
not material coefficients or point containment. The feature remains opt-in;
default behaviour and production runtimes are unchanged.

## Spatial queries without full preparation

Requires CSXCAD with the public `CSPrimitiveIndex(structure, config)` constructor.
This is a C++ setup acceleration experiment, not a new geometry, XML, or physics
option.
It is off by default. Both scripting interfaces inherit the environment of
the launched solver:

```bash
OPENEMS_GEOMETRY_INDEX=16,16,1 openEMS model.xml --numThreads=4
```

In Python, set `os.environ['OPENEMS_GEOMETRY_INDEX'] = '16,16,1'` before
`FDTD.Run(...)`; in Matlab/Octave, use
`setenv('OPENEMS_GEOMETRY_INDEX', '16,16,1')` before `RunOpenEMS(...)`.
Unset the variable or use `0` for the legacy path. Invalid configuration and
requested support missing from CSXCAD are errors, not silent fallbacks.

The three integers select logical X/Y/Z tiles in the mesh domain. They do
not refine or change the Yee grid. `Nz=1` suits the existing full-Z column
queries. More tiles can reduce candidates but increase construction time and
storage; measure on representative inputs. Each of the three type-specific
indices is limited to 1,048,576 buckets and 8,388,608 stored rank references.
These are storage-count limits, not a process RSS guarantee; vector overhead,
capacity and transient query results consume additional memory.

CSXCAD builds the immutable, non-owning index and answers candidate queries;
openEMS chooses the bounds/bins and lifecycle. Separate material,
material-plus-metal and excitation indices retain the legacy type-specific priority/tie
order. Exact material containment and PEC decisions are unchanged. Transforms
and unknown bounds remain conservative candidates. Wide/outside queries reuse
the sorted snapshot. Other property masks retain the legacy query path.

The indices cover material/PEC coefficient construction, extension construction
(including UPML, dispersive materials and excitation) and setup debug dumps.
They are released before time stepping, including early
return and exception paths. The geometry must have completed `Update()` and
must remain unchanged while building the operator, as for the parallel legacy
setup. The built-in extension builders do not mutate CSXCAD geometry. Derived
operators/extensions that mutate geometry during any of those stages
are incompatible with this opt-in path and must leave it disabled. Virtual
`GetPrimitivesBoundBox` dispatch is retained.

Rebuild openEMS when changing the CSXCAD headers/library; do not mix a new
header with an older shared library. This branch supports building against an
older CSXCAD too, with indexing disabled. It does not qualify an existing
production binary or experimental material extension automatically.

With the companion polygon-bounds revision in CSXCAD, untransformed Cartesian
`Polygon` and `LinPoly` figures are spatially culled using conservative bounds
expanded by drawing tolerance. Unsupported/transformed cases remain conservative.
Legacy candidate lists may contain additional distant polygons; exact material
evaluation, priorities and coefficients must remain equivalent. Rebuild native
dependents with the matching CSXCAD headers and library.

CSXCAD also offers optional borrowed, precomputed tile lists. Spatial-only mode
uses filtered `Query()` results: broader tile lists can cost more during the
much larger number of exact point tests. Prepared mode stores those exact
filtered results in advance; it does not substitute broader tile candidates.
No per-query shared timer, mutable cache or FDTD-loop change is introduced.

## Public local benchmark (14 November 2026)

Unchanged upstream `benchmark/bench.py`, public companion `geometry_index.py`,
three interleaved master/prepared pairs, XS/basic/two threads/0.5 s target.
Both binaries use the same new CSXCAD. The synthetic polygon case adds 1024
extruded polygons; no private geometry is included.

| Case | Master setup | Prepared setup | FDTD master/prepared MC/s | Spread | Probes |
|---|---:|---:|---:|---:|---|
| Inhomogeneous | 0.03551 s | 0.04238 s | 267.3 / 267.2 | 0.5% | identical |
| Polygons | 0.58375 s | 0.08361 s | 269.4 / 267.6 | 1.2% | identical |

Setup is the benchmark's wall-minus-FDTD metric. The polygon comparison includes
conservative polygon culling as well as prepared queries; it does not isolate
the table change. The tiny control pays 6.86 ms extra. FDTD differences are
classified as noise. This is local CPU evidence, not full-antenna or device
qualification. Release/ASan/UBSan and old-CSXCAD controls pass; indexed/default
coefficient arrays and 48-step basic-engine fields are identical.
