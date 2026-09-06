# Seeding Strategy

Seeding is **centralized in `runTractography.m`**. It selects the seed mask, and the
trackers only consume it: `nim_tractography_standard` and `nim_tractography_hinec` validate
that `options.seed_mask` is present and error out if it is not, rather than falling back to a
mask of their own. Policy (what to seed) and mechanism (how to track) stay separated, which
is what lets challenge-specific code and ROI experiments change the seeds without touching a
tracker.

```
runTractography.m                         nim_tractography_*.m
     │                                             │
     ├─ load the nim                               │
     ├─ pick a seeding strategy                    │
     ├─ build seed_mask, set options.seed_mask ────┼─► validate seed_mask is non-empty
     │                                             ├─► place sub-voxel seeds (nim_seed_offsets)
     └─ dispatch on algorithm ────────────────────►└─► track bidirectionally from each seed
```

---

## Strategy hierarchy

`runTractography` tries these in order and reports which one it used.

### 0. Explicit ROI — `seeding.roi`

Highest priority. Seeds go only in the named atlas regions, resolved by
`nim_roi_mask` from JHU indices and/or names (freely mixed), optionally dilated by
`seeding.roi_dilate` voxels.

```yaml
tractography:
  seeding:
    roi: [SLF_L, SLF_R]
    roi_dilate: 1
```

The ROI is then intersected with the brain mask and the FA floor exactly as the other
strategies are, so switching to ROI seeding changes *where* seeds go and nothing else. An
ROI that survives none of those intersections is an error, reported with the voxel count
after each stage so it is clear which one emptied it.

This is the cheapest way to interrogate a single bundle, and it is what bundle-specific
reconstructions and fixed-seed convergence ladders use.

### 1. Preprocessed brain mask

Used when `nim.mask` is present and non-empty: `nim.mask > 0.5`, refined by the FA floor.
Whole-brain coverage with clean boundaries — the normal path after `main.m`.

### 2. Expanded parcellation mask

Used when there is no brain mask but `nim.parcellation_mask` exists: the labelled regions
dilated by a 3-voxel spherical structuring element to pick up the surrounding white matter,
then constrained to `FA > 0.05`.

### 3. FA-threshold fallback

Used only when neither mask exists: `FA > 0.10`. The run warns, because this covers the
white-matter core and misses low-anisotropy structures such as the fornix and cingulum.

### The FA floor

`seeding.fa_min` (default `0.05`) sets the minimum FA for a voxel to be seeded. The default
excludes CSF and little else. Raising it to ~0.2 gives white-matter-only seeding, which is
what a matched WM-seeding comparison needs.

---

## Seed placement within a voxel

`seeding.density` (default `8`) seeds per seeded voxel, placed by `nim_seed_offsets` on a
deterministic sub-voxel lattice in \([-0.5, +0.5]^3\):

- a **perfect cube** (1, 8, 27, 64, …) gives the full `per_axis³` lattice;
- any other count gives a deterministic farthest-point subset of the next larger lattice, so
  the seeds stay spread through the voxel instead of clumping in a corner.

The count is honoured **exactly**. An earlier implementation rounded up to the next perfect
cube, so `density: 4` silently placed 8 seeds per voxel. There is no RNG in the `uniform`
strategy, so runs are reproducible and streamlines can be compared one-to-one across a
parameter ladder. `seeding.strategy: random` jitters the placement instead.

Tracking is bidirectional from each seed, so the upper bound on tracks is
`2 × density × seed voxels`, before any track is discarded by `termination.min_arc` or the
ROI filters.

---

## Where a track may go is not where it may start

The **propagation mask** is derived independently of the seed mask, in both the `hinec` and
`standard` trackers:

```
prop_mask = (FA > termination.fa_min) & (brain mask)
```

dilated by one voxel for the boundary test, and overridable through
`options.propagation_mask`. It used to be `imdilate(seed_mask, ...)`, which confined every
streamline to a one-voxel skin around the seed region. With whole-brain seeding the two
coincide and nothing looks wrong; with ROI seeding the propagation domain collapsed onto the
ROI and tracks died almost immediately. Keep the two distinct when modifying either tracker.

---

## Special case: IronTract

The IronTract Challenge seeds from an injection site that may lie outside the brain mask.
This is handled downstream in `src/nim_challenges/nim_irontract_submit.m`, which unions the
injection mask into the seed mask, rather than by adding a strategy to `runTractography`.

---

## Checking what a run actually seeded

The seeding block prints the chosen strategy, the seed-voxel count and its share of the
volume, the seeds-per-voxel figure actually realised, and the approximate inter-seed
spacing. To compare strategies, remove the brain mask from a copy of the nim and re-run:

```matlab
S = load('data.mat');
S.nim.mask = [];              % force the FA-threshold fallback
save('data_no_mask.mat', '-struct', 'S');
runTractography('data_no_mask.mat');
```

A brain-mask seeding run and an FA-threshold run differ in both the number and the
distribution of seed voxels; the FA floor is the reason low-anisotropy bundles drop out of
the fallback. Compare the two the same way as any other experiment — one scorable run
directory each, compared by `scoring/renauld2023/results.json`. See
[Tractography Methods](TRACTOGRAPHY_METHODS.md#config-driven-experiment-workflow).

---

## Common issues

| Symptom | Cause | Fix |
|---|---|---|
| `No seed mask provided` | a tracker was called directly, without `options.seed_mask` | go through `runTractography` |
| `ROI seed mask is empty after masking` | the ROI vanished under the brain mask or the FA floor | raise `seeding.roi_dilate`, lower `seeding.fa_min`, or check the region names |
| Low-anisotropy bundles missing | the FA-threshold fallback was used | make sure `nim.mask` is present in the nim |
| Fewer or more tracks than expected | `seeding.density` is per **seeded voxel**, and tracking is bidirectional | check the printed seed-voxel count |

---

## See also

- [Tractography Methods](TRACTOGRAPHY_METHODS.md) — the trackers and the experiment workflow.
- [YAML Config](YAML_CONFIG.md) — `seeding.*` and `filter.*` in full.
- [IronTract Workflow](IRONTRACT_WORKFLOW.md) — injection-site seeding end to end.
