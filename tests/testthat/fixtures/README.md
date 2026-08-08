# PhysioAnalysis test fixtures

## `idw_reference.rds` — inverse-distance-weighting topomap grid (WSCB-06, WSCB-09)

A closed-form interpolation reference used to cross-validate the topographic-map
interpolation `.interpolateIDW()` (the correction in WSCB-06 documents it as
inverse-distance weighting / Shepard's method, not Perrin spherical splines).

- **Reference method**: the Shepard inverse-distance-weighting estimate
  (power 2) on a synthetic 9-electrode montage (an eight-point ring plus the
  centre of the unit disc), computed by an **independent** vectorised route
  (a full weight matrix rather than the package's per-point loop) with the
  exact-at-nodes rule applied explicitly.
- **Contents**: `x`, `y`, `values` (electrode positions and values), `power`,
  `grid_x`, `grid_y` (a 21×21 evaluation grid), `expected` (the reference
  interpolated values), `provenance`.
- **Gate** (`test-topomap.R`): `.interpolateIDW()` reproduces `expected` to
  within 1e-10 across the grid, and returns each electrode's own value exactly
  when evaluated at the electrode positions (the Shepard interpolation
  property).

Regenerate with:

```r
Rscript physio-ecosystem/publication/scripts/validate_correctness_fixes.R --write
```

## `spherical_spline_mne_reference.rds` - Perrin interpolation (WSCB-10)

An offline reference generated with MNE-Python's
`mne.channels.interpolation._make_interpolation_matrix()` from a pinned MNE
release. The fixture uses the same shared upper-hemisphere lift as
PhysioAnalysis before calling MNE, with `stiffness = 4`, 50 Legendre terms, and
`alpha = NULL`.

- **Contents:** planar source/query coordinates, shared center/radius, source
  values, MNE interpolation matrix output, parameters, and provenance.
- **Gate:** normalized RMS error below 0.05, with a tighter implementation
  comparison recorded by the test.
- **Generation:** Python and MNE version, source URL, exact command, and
  generation date are stored in the fixture's `provenance` list.

The fixture contains numeric inputs and output only. MNE is not a package or
test dependency.

Regenerate from the `PhysioAnalysis` package directory after installing the
pinned MNE wheel into `/tmp/mne-wscb10`:

```bash
_MNE_FAKE_HOME_DIR=/tmp/mne-home \
  PYTHONPATH=/tmp/mne-wscb10 \
  python3 tests/testthat/fixtures/generate_spherical_spline_mne_reference.py \
  /tmp/wscb10_mne_fixture.json
Rscript -e 'x <- jsonlite::read_json("/tmp/wscb10_mne_fixture.json", simplifyVector=TRUE); saveRDS(x, "tests/testthat/fixtures/spherical_spline_mne_reference.rds", version=3)'
```
