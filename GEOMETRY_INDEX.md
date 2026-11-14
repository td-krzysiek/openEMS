# Experimental setup geometry index

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

CSXCAD also offers optional borrowed, precomputed tile lists. This consumer
uses filtered `Query()` results: broader lists can cost more during the much
larger number of exact point tests. No per-query shared timer, mutable cache or
FDTD-loop change is introduced.
