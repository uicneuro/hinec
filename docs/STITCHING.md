# Tractlet stitching (experimental)

`algorithm: stitching` generates dense, short, bidirectional tractlets and joins
compatible endpoints into longer chains. The implementation is
`nim_tractography_stitching`.

## Design and implementation plan

1. Generate the requested deterministic subvoxel seed lattice. Start one tractlet
   per seed for DTI; for CSD, start one for each available seed-voxel peak.
2. Track only `stitching.fragment_arc` total voxel arc, split between two arms.
   Direction mode uses midpoint RK2 on local, trilinearly interpolated dyadics.
   CSD selects the nearest incoming-direction peak at each of the eight corners
   before interpolating its dyadic. There is no long-range spline stencil.
3. With `stitching.geometry: mmf`, evolve a fresh moving frame per tractlet with
   midpoint integration of the supplied curvature/torsion connection. DTI reads
   both curvature and torsion. CSD reads per-peak curvature; the current geometry
   provider has no per-peak torsion, so no torsional twist or torsion join gate is
   applied in that mode. Curvature clipping is explicit and configurable.
4. Form a bounded nearest-neighbor endpoint graph. A join needs a short gap,
   compatible tangents, forward extension from both endpoints, and a sampled
   cubic Hermite bridge supported by the propagation mask, FA, and field direction.
   MMF also checks endpoint curvature and, where defined, torsion consistency.
5. Process edges by increasing gap/alignment/geometry cost. Each endpoint is
   consumed at most once. Union-find rejects graph cycles; component arc limits
   prevent unbounded chains. This is deterministic greedy assembly, not a global
   optimum or a branching graph search.
6. Assemble the chains and apply `termination.min_arc` to their actual polyline
   lengths. Validate synthetic joins/crossings/gaps, then score real ROI pilots.

Distances and curvature units follow the existing tracker contract: voxels and
inverse voxels. At 2 mm isotropic resolution the default 3-voxel fragment is
6 mm long, and the 1-voxel join radius is 2 mm. `termination.max_arc` caps a
whole assembled chain for this algorithm. `integrator.step` controls integration
and bridge sampling. Supported settings are `rk2`, `trilinear`, `upsample: 1`,
`act: false`; unsupported settings fail explicitly. MMF anchoring is not used:
every short fragment starts with a fresh field-aligned tangent and reference frame.

## Run

```bash
./bin/run_tractography.sh stitching_dti --score --set seeding.roi=Fornix
./bin/run_tractography.sh stitching_csd --score --set seeding.roi=Fornix
./bin/run_tractography.sh stitching_mmf_dti --score --set seeding.roi=Fornix
./bin/run_tractography.sh stitching_mmf_csd --score --set seeding.roi=Fornix
```

The presets use eight seeds per voxel. They are intended for initial ROI
experiments: the explicit `stitching.max_fragments` budget fails before generation
if seed count times possible peaks exceeds the budget. Increase it deliberately
or reduce density/ROI. Endpoint search is bounded by `stitching.neighbors`; a
small value can miss compatible candidates in very dense neighborhoods.

A fragment-only ablation is `--set stitching.join=false`. With the normal
minimum final length this may correctly produce zero tracks (the pipeline then
reports no tracks). Use `--set termination.min_arc=0` to save and inspect those
short fragments; scores with that changed length threshold are not a controlled
comparison of final long-tract quality.

## Inspect and evaluate

`track_meta.graph` contains candidate and accepted edge counts, bridge rejection
counts, discarded short-chain count, accepted endpoint IDs and costs, and the
source fragment IDs of every retained chain. `fragment_seed_indices` maps each
retained chain to all its source seeds. `seed_index`/`seed_points` identify its
first source fragment, not a claim that a chain originated from only one seed.
With `debug.trace: true`, all generated fragments, their seeds, and MMF endpoint
geometry are also saved. This trace can be large; `debug.trace_max` does not cap
stitching fragment provenance.

Use synthetic straight, curved, crossing, and masked-gap cases before relying
on bundle scores. Compare matched ROI runs against HINEC/MMF and GT for geometry,
coverage, invalid connections, and runtime. High valid counts alone do not prove
better anatomical reconstruction. Stitching can connect the wrong nearby fibers,
and short integration cannot recover information absent from the estimated field.
Graph cycle prevention does not guarantee that a geometric chain has no spatial
self-intersections. No GT fibers or scorer gates are consulted while joining.

## Dense-seeding experiment controls

`--set stitching.search=forward` searches all endpoints in the forward cone
within `join_radius`, applies tangent and MMF geometry gates, and then limits
compatible candidates to `neighbors`. The limit is distributed round-robin
across eight strata: four transverse azimuth quadrants and two gap-distance
shells. Within each stratum candidates are ranked by the usual join cost, with
endpoint ID breaking ties. The final greedy assembly and bridge checks are
unchanged. The default `knn` search remains available for controlled comparisons.

Radius queries are batched (256 endpoints), and candidate storage is bounded by
`2 * fragments * neighbors` before pair deduplication. Computation still grows
with the number of endpoints inside the query sphere; this is not a constant-cost
search. The enclosing sphere contains the full allowed forward cone; compatibility
checks remove the extra sphere volume. Seed grids of 8, 64, and 216 use respectively
2, 4, and 6 evenly spaced subdivisions on each voxel axis. Set `max_fragments`
to accommodate the seedable voxel count multiplied by the density. CSD budgets
also include the maximum possible seed-voxel peak count.

Compare densities within each algorithm, and compare against continuous HINEC/MMF
at the same densities. Seed density and endpoint search strategy are separate
experimental factors. ROI/whole-brain identity for independent seed trackers does
not extend to stitching: adding outside fragments can change the greedy graph.
