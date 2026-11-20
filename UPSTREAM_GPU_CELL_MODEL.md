# GPU domain integration with upstream master

The source includes upstream master `c22a063841d6f09d13504670ae119b4dca702cdb`
(19 November 2026). The integration combines that tree with the local-domain
runtime and setup improvements from the fork. It does not replace an installed
or running binary. Use matching core, GPU plugins and rebuilt Python bindings.

## Changes from upstream

The upstream GPU engine now expands the CPU cell model instead of uploading
four dense coefficient arrays. UPML and supported dispersive materials belong
to the main update. This integration partitions that expanded model into owned
X slabs with a field halo on each side. Each slab has compact combination IDs
and history slots for owned nodes only. Mesh factors and material profiles
retain the global discretization.

Point operations keep global coordinates for hook-box filtering and remap their
field indices to local storage. Excitation, lumped RLC, remaining dispersive
orders and experimental opaque surfaces execute only on their owning domain.
The opaque-surface operation uses the common reference/HIP/Metal point kernel
interface. Models with more poles than the cell model supports (including the
eight-branch LossyMetal model) retain the generic device-side ADE hooks.

Halo exchange runs after the native correction/source hooks. Range readback,
owned-only energy reductions, explicit cuts, sampled domain profiling and HIP
serial/async/overlap transport remain part of the domain adapter. Energy sums
accept global prefix bounds and exclude both halo planes. Unsupported domain
extensions fail explicitly.

The fork's surface models now use upstream `GetEC`/`SetEC` accessors. Direct
writes to the retired dense circuit arrays are invalid when upstream's distinct
material table is active. The numerical models and stabilization rules are
unchanged by this API migration.

## Explicit limits

- Multi-domain execution uses separate E/H steps. Upstream temporal blocking
  and fused E/H trial selection remain available for a single executor; they
  are disabled across domains until halo-aware schedules are qualified.
- `OPENEMS_GPU_MAIN_COEFF=auto` selects the cell model. The former `raw` mode
  fails explicitly. `OPENEMS_GPU_LOSSY_FUSION` must be unset; its former branch
  kernel A/B controls are not implemented by the cell-model runtime. The
  archived policy helpers/tests do not qualify the new update path.
- Cartesian multi-domain support does not establish cylindrical or Cartesian
  subgrid support. This change does not extend the separate subgrid prototype.
- Logical Metal domains exercise ownership and queue ordering on one physical
  GPU; they do not establish physical multi-GPU scaling or transport behavior.
- A successful small-fixture comparison is not a speed, large-model memory,
  numerical convergence or application-validation claim.

## Reproducible checks

Configure with `ENABLE_DOMAIN_RUNTIME=ON`, `OPENEMS_DOMAIN_TESTS=ON` and
`OPENEMS_BUILD_TESTS=ON`, then build and run `ctest --output-on-failure`.
The native reference copy and energy tests run without a GPU or Python binding.

```sh
python3 tests/manual/gpu_cell_domains.py --binary /candidate/bin/openEMS \
  --output /new/reference-cell-domains
python3 tests/manual/gpu_cell_domains.py --binary /candidate/bin/openEMS \
  --engine gpu --halos serial --output /new/physical-cell-domains
```

The fixed-step matrix covers ten public fixtures and four layouts: one
executor, unequal two-domain cuts, a one-line middle domain through the source
and conductor, and one-line domains at both physical boundaries. It compares
all complex E/H dump components and the cut-plane probe against the same
build's single-executor result, with a peak-relative tolerance of `1e-6`.
Run `--halos async` and `--halos overlap` on physical HIP devices separately.

The existing `tests/manual/domain_backend.py` also covers equal 2/3/8-domain
layouts; it requires matching CSXCAD Python bindings with LossyMetal support.
The older coefficient/fusion A/B harnesses target the previous runtime and
must not be used as evidence for this integration.

## Recorded qualification

The 20 November qualification covers a native Metal/domain build and a Linux
HIP-as-CUDA plugin build for architectures 86/89, each with 14/14 CTest checks.
Fresh Python bindings pass all 19 unit-test modules. The fixed-step matrices
cover 40 reference runs and 40 Metal runs (ten fixtures, four layouts); compared
complex fields and probes are identical within each backend. The material-table
fallback and two-domain checkpoint/snapshot continuation controls also pass.
Legacy control rejection is checked against the native reference engine.

These are local/public-fixture controls. Physical multi-device HIP transport
and performance remain unqualified: the physical devices were leased by other
work during this integration. The qualification does not authorize replacing
an existing runtime. Offline Metal compilation needs the optional Apple Metal
Toolchain; runtime compilation and execution are tested on the actual device.
