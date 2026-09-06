# ISMRM 2015 Challenge Scoring

The ISMRM 2015 Tractography Challenge ships a synthetic whole-brain phantom
together with the ground truth used to score it. This page documents what is in
that scoring package, how a bundle is defined by it, what the shipped HINEC
scoring path does with it, and what our reconstructions currently measure.
Everything below is read from
`data/ismrm2015/scoring_data_Renauld2023/` (Renauld et al. 2023 revision of the
original 2015 scoring system).

!!! warning "The scorer expects a whole-brain tractogram"
    The Renauld 2023 scorer takes **one whole-brain tractogram** and segments it
    into 26 bundles using ROI gates. Its headline numbers — `mean_f1`,
    `mean_OL`, `mean_OR_gt`, `VB` — are averages over all 26. Feeding it
    streamlines from a single seed ROI is not a valid use of it: at most one
    bundle can be recovered, so `VB` caps at 1 of 26 and a mean over 26 bundles
    is dominated by the 25 that were never attempted. **No score obtained that
    way is comparable to a published challenge result.** Every score collected
    for HINEC so far is of that kind; the whole-brain submission has not been
    run.

---

## Reconstruction against ground truth

Each bundle below is seeded from its own ISMRM mask, tracked with **no
filtering**, and then **clipped** to the bundle's containment corridor for
display. Grey is the ground-truth tractogram; colour is ours.

!!! warning "Clipped, not filtered — and why that distinction matters"
    The scorer defines a bundle by containment: a streamline with *any* point
    outside the corridor is not that bundle and is discarded whole. That is the
    correct rule for scoring and a misleading one for a picture. A streamline
    that follows a bundle correctly and then overshoots past its end gets deleted
    entirely, taking the correct portion with it — so the figure shows a bundle
    that stops short, which reads as the tracker *failing to reach* when what
    actually happened is that it reached too far.

    These figures clip instead, and state how much length was lost, so the
    overshoot is visible rather than hidden. Clipping is `nim_plot_bundles`'
    `.clip` option; rejection is `nim_filter_tracks_roi`, which is what scoring
    uses.

**How much of each bundle leaves its own corridor:**

| bundle | streamlines | length outside the corridor |
|---|--:|--:|
| CC_u_shaped | 16287 | **28%** |
| ILF_right | 2394 | **42%** |
| BPS_right | 7810 | **43%** |
| Cingulum_right | 5450 | **59%** |
| SLF_right | 3901 | **61%** |
| UF_right | 1398 | **62%** |

That is the honest headline for single-tensor DTI tracking on this phantom:
between a quarter and two thirds of produced streamline length leaves the bundle
it was seeded in. The streamlines are typically correct along the bundle and then
continue past its end onto whichever tract is locally strongest — a crossing
problem, not a termination one. `field: csd` is the intended remedy and has not
yet been evaluated here.

These remain **qualitative** figures, and seeding from a bundle's own mask is a
far easier problem than the whole-brain submission the challenge is designed
around, so none of this constitutes a challenge score.

### Projection

![BPS_right](img/bundle_bps_right.png)

### Commissural

![CC_u_shaped](img/bundle_cc_u_shaped.png)

### Association

![SLF_right](img/bundle_slf_right.png)

![ILF_right](img/bundle_ilf_right.png)

![UF_right](img/bundle_uf_right.png)

### Limbic

![Cingulum_right](img/bundle_cingulum_right.png)

---

## What defines a bundle

A bundle in this scoring system is an **endpoint pair plus a containment
corridor**, not a parcellation label. `config_file_segmentation.json` gives one
entry per bundle:

```json
{
  "UF_right": {
    "all_mask": "ROI/all_masks/UF_right.nii.gz",
    "head":     "ROI/endpoints/UF_right_head.nii.gz",
    "tail":     "ROI/endpoints/UF_right_tail.nii.gz"
  }
}
```

A streamline is assigned to `UF_right` only if one end lands in the head ROI, the
other end in the tail ROI, and **every** point lies inside `all_mask`. Two
optional gates appear on some bundles:

- `any_mask` — at least one point must fall inside it. Present on 6 of the 26
  entries: `CC_u_shaped`, `MCP`, and the four `ICP_*_part*` entries.
- `length`, `length_x`, `length_y`, `length_z`, `length_x_abs` — total streamline
  length and per-axis extents, in millimetres. Present on 4 of the 26:
  `CC_u_shaped`, `Cingulum_left`, `Cingulum_right`, `MCP`. `CC_u_shaped` is the
  fully loaded case:

```json
{
  "CC_u_shaped": {
    "all_mask":     "ROI/all_masks/CC_u_shaped.nii.gz",
    "any_mask":     "ROI/any_masks/CC_u_shaped_inclusion_mask.nii.gz",
    "head":         "ROI/endpoints/CC_striatal_left.nii.gz",
    "tail":         "ROI/endpoints/CC_striatal_right.nii.gz",
    "length":       [70, 1000],
    "length_y":     [0, 32],
    "length_x_abs": [35, 1000]
  }
}
```

!!! warning "An atlas label of the same name is not the same region"
    A JHU atlas label and an ISMRM bundle can share a name and describe very
    different volumes. JHU label 47, *Uncinate fasciculus R*, is 24 voxels on
    the 2 mm DWI grid; the corresponding ISMRM bundle-density mask (`scoring_data_2015/masks/bundles/`) occupies 1503 and the Renauld containment corridor 14260, and the two
    overlap at a Dice of 0.018. Seeding from the atlas label and scoring against
    the ISMRM bundle of that name is therefore not a like-for-like comparison.
    Address the challenge's own regions instead: build the parcellation from the
    challenge masks with `nim_parcellation_from_masks`, after which ROI names
    resolve against `nim.roi_masks` (see `nim_roi_mask`).

### Bundle gates in HINEC

Two `tractography.filter` predicates express the scorer's definition directly:

| key | test |
|---|---|
| `filter.endpoints_in` | two regions; keep a track only if one **end** lands in each, either way round |
| `filter.contained_in` | keep a track only if **every** point lies inside the given regions |

These are distinct from `filter.include_roi`, which is a *waypoint* test — it
asks whether a track passes through a region, not where it stops.

```yaml
tractography:
  filter:
    endpoints_in: [UF_right_head, UF_right_tail]
    contained_in: [UF_right]
```

---

## What is in the scoring package

```
scoring_data_Renauld2023/
├── bundles/                       21 ground-truth .trk files
│   └── sub_bundles/               8 further .trk files (CC and ICP subdivisions)
├── ROI/
│   ├── all_masks/                 26 containment corridors
│   ├── any_masks/                 4 inclusion masks
│   └── endpoints/                 45 head/tail ROIs (some shared, e.g. brainstem)
├── config_file_segmentation.json  26 bundle definitions (ROI gates)
├── config_file_tractometry.json   bundle name -> ground-truth .trk
└── t1.nii.gz                      1 mm reference, 180x216x180
```

The 21 top-level bundles are commissural (`CA`, `CC`, `CP`, `MCP`), association
(`Cingulum`, `ILF`, `OR`, `SLF`, `UF`, each left and right), projection (`BPS`,
`ICP`, `SCP`, each left and right) and `Fornix`. The segmentation config reaches
26 entries because `CC` is scored as four sub-bundles (`CC_temporal`,
`CC_u_shaped`, `CC_ventro_striatal1`, `CC_ventro_striatal2`) and each `ICP` as
two parts.

Ground-truth streamline counts span two orders of magnitude, which is why a mean
over bundles is not a mean over streamlines:

| bundle | streamlines | | bundle | streamlines |
|---|--:|---|---|--:|
| MCP | 21008 | | ILF_left | 11164 |
| Cingulum_right | 20618 | | BPS_left | 11162 |
| CC | 16550 | | ILF_right | 10630 |
| BPS_right | 15400 | | OR_right | 9524 |
| Cingulum_left | 14206 | | OR_left | 7252 |
| SLF_left | 12483 | | UF_right | 6751 |
| SLF_right | 11920 | | UF_left | 5899 |
| ICP_left | 4217 | | Fornix | 3827 |
| ICP_right | 3224 | | SCP_left | 1795 |
| SCP_right | 1560 | | CA | 430 |
| CP | 365 | | | |

`CP` and `CA` are the smallest and, together with the commissures, the hardest to
recover.

---

## The scoring path

`bin/run_ismrm_scoring.sh <run_dir>` is the single entry point. It performs three
steps:

1. **Convert.** `scripts/hinec_to_trk.py` turns the newest
   `<run_dir>/tractography/tracks*.mat` into `scoring/tracks.trk` in RAS world
   millimetres, using the preprocessed DWI affine from `<run_dir>/intermediate/`
   for voxel→world placement (falling back to `data/ismrm2015/ismrm2015.nii.gz`)
   and attaching the scoring `t1.nii.gz` as the saved TRK reference so scilpy's
   space-compatibility check against the ROI masks passes.
2. **Score (Renauld 2023).** `scil_tractogram_segment_with_ROI_and_score` runs
   against `config_file_merged.json` — the segmentation rules and the per-bundle
   `gt_mask` in one file, produced by
   `scripts/build_ismrm_scoring_config.py`. Without the merged config the run
   yields bundle counts but no Dice or overlap. Output:
   `<run_dir>/scoring/renauld2023/results.json`.
3. **Cross-check (optional).** If `data/ismrm2015/scoring_data_2015/` and the
   original challenge scorer are present, the dedicated 2015 scorer is run as
   well, into `<run_dir>/scoring/dedicated2015/`.

Headline keys in `results.json`:

| key | meaning |
|---|---|
| `mean_f1` | Dice-style agreement per bundle, averaged over the 26 |
| `mean_OL` | overlap — fraction of the true bundle covered (recall) |
| `mean_OR_gt` | overreach — produced volume outside the true bundle |
| `VB` | valid bundles: how many of the 26 were recovered at all |

Bundle recognition by shape (RecoBundles) is an alternative the challenge
tooling supports and the HINEC path does not use; it is more forgiving of ROI
misalignment and a correspondingly weaker test of anatomical placement.

---

## Coordinate spaces

Scores collapse to zero on a space mismatch, and a space mismatch looks exactly
like a tracking failure, so it is worth checking before believing a low score.

| | HINEC tracks | ISMRM ground truth |
|---|---|---|
| format | `.mat`, cell array of N×3 | `.trk` (TrackVis) |
| coordinates | voxel indices, 1-based (MATLAB) | RAS world millimetres |
| grid | DWI, 2 mm, 90×108×90 | T1, 1 mm, 180×216×180 |
| space attribute | implicit | explicit in the TRK header |

Since August 2022 the challenge scoring code loads tractograms through Dipy's
`StatefulTractogram` and no longer applies the half-voxel shift that older
TrackVis-era readers used. Coordinates are taken from the header as given, so a
TRK written with the wrong reference is silently misinterpreted rather than
rejected: the ROI gates then match nothing and every bundle scores 0.

`hinec_to_trk.py` handles the conversion, and `nim_read_trk` performs the reverse
mapping (world RAS → DWI voxel space) when ground truth is loaded into MATLAB for
plotting. Skipping that reverse step draws the ground truth at half scale in a
corner — which, again, reads as a tracking failure and is a units bug.

To check alignment directly, compare coordinate ranges and header spaces of the
two tractograms:

```python
import nibabel as nib
import numpy as np

ours = nib.streamlines.load('tracks.trk')
gt   = nib.streamlines.load('scoring_data_Renauld2023/bundles/CA.trk')

for name, trk in [('ours', ours), ('ground truth', gt)]:
    pts = np.vstack(list(trk.streamlines))
    print(name, [f'[{pts[:, i].min():.1f}, {pts[:, i].max():.1f}]' for i in range(3)])
```

Both should span the same millimetre range. If they do not, fix the reference
before reading anything into the scores.

---

## Limits of this benchmark

- **One phantom, one subject.** Nothing measured here establishes generality.
- **The ground truth is synthetic.** The bundles are anatomically motivated but
  simulated; agreement with them is not agreement with a real brain.
- **ROI-gate sensitivity.** Bundle assignment depends on precise mask alignment,
  so registration error is charged to the tracker.
- **Bundle-averaged metrics.** `mean_f1` weights `CP` (365 ground-truth
  streamlines) as heavily as `MCP` (21008).

For a complementary check with a biological rather than synthetic reference, see
[IronTract](IRONTRACT_WORKFLOW.md), which scores against tracer injections.

---

## References

- [ISMRM 2015 Tractography Challenge](http://tractometer.dinf.usherbrooke.ca/ismrm_2015_challenge/)
- [Renauld 2023 scoring code](https://github.com/scilus/ismrm_2015_tractography_challenge_scoring)
- [Dipy StatefulTractogram](https://dipy.org/documentation/latest/reference/dipy.io.stateful_tractogram/)
- [TrackVis file format](http://trackvis.org/docs/?subsect=fileformat)
