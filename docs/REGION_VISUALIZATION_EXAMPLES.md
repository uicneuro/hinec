# Region-Specific Tractography Visualization

Worked examples for displaying the streamlines associated with one or more
parcellation regions. All of it goes through `visualizeTractography` in
`mode: 'region'` — the older standalone entry points
(`visualizeTractographyRegion`, `visualizeTractographyAllRegions`,
`nim_plot_tractography_region`) were folded into it and no longer exist as
callable functions.

For the full parameter reference, see the
[Visualization Guide](VISUALIZATION_GUIDE.md).

---

## Prerequisites

You need a processed `nim` with a parcellation and a set of tracks:

```bash
./bin/run_hinec.sh data/parcellation_sample/sample sample.mat config/hinec_default.yml
./bin/run_tractography.sh hinec_dti     # defaults to the newest processed run
```

Both land in a timestamped `hinec_runs/run_<TIMESTAMP>_<config>/`, with the nim
under `output/` and the tracks under `tractography/` as
`tracks_<algorithm>_<timestamp>.mat`.

---

## Basic usage

`visualizeTractography` accepts either a run directory (it finds the nim and the
newest tracks itself) or an explicit tracks/nim pair.

```matlab
% From a run directory
visualizeTractography('hinec_runs/run_<TIMESTAMP>_hinec_dti/', ...
    'mode', 'region', 'region', 5);

% From explicit files
visualizeTractography('hinec_runs/run_<TIMESTAMP>_hinec_dti/tractography/tracks_hinec_<TS>.mat', ...
                      'hinec_runs/run_<TIMESTAMP>_hinec_dti/output/sample.mat', ...
    'mode', 'region', 'region', 5);
```

The parameter is `region` (singular), and it takes a scalar or a vector:

```matlab
visualizeTractography(run_dir, 'mode', 'region', 'region', [5, 10, 15]);
```

---

## Filter modes

`filter_mode` decides what "associated with the region" means:

| mode | keeps tracks that |
|---|---|
| `'pass_through'` (default) | touch the region at any point |
| `'start_in'` | begin inside the region |
| `'end_in'` | end inside the region |
| `'connect_to'` | link the region to another parcellation region |

```matlab
visualizeTractography(run_dir, 'mode', 'region', 'region', 12, 'filter_mode', 'start_in');
visualizeTractography(run_dir, 'mode', 'region', 'region',  8, 'filter_mode', 'end_in');
visualizeTractography(run_dir, 'mode', 'region', 'region',  5, 'filter_mode', 'connect_to');
```

Two thresholds shape `pass_through`: `min_overlap` (fraction of the track inside
the region, default 0.05) and `min_region_points` (absolute number of in-region
points, default 3). Raise either to demand more than a glancing intersection.

!!! note "Display filtering is not tractography filtering"
    These options select what is *drawn* from a track set that already exists.
    To restrict which streamlines are *produced*, use the tractography config —
    `seeding.roi` to seed inside a region, `filter.include_roi` for a waypoint
    test, and `filter.endpoints_in` / `filter.contained_in` for the endpoint and
    containment tests the ISMRM scorer uses. See
    [ISMRM Scoring](ISMRM_SCORING_ANALYSIS.md).

---

## Colouring

| `color_mode` | effect |
|---|---|
| `'direction'` (default) | RGB from local orientation: red = left–right, green = anterior–posterior, blue = superior–inferior |
| `'fa'` | hot colormap from mean FA along the track |
| `'uniform'` | one colour for all tracks; the cheapest to render |
| `'region'` | one colour per parcellation region |

```matlab
visualizeTractography(run_dir, 'mode', 'region', 'region', 5, 'color_mode', 'fa');
```

---

## A fuller example

```matlab
visualizeTractography(run_dir, ...
    'mode', 'region', ...
    'region', 5, ...
    'filter_mode', 'pass_through', ...
    'color_mode', 'fa', ...        % colour by mean FA
    'max_tracts', 500, ...         % cap the number drawn
    'show_region', true, ...       % draw the region as a 3D overlay
    'region_alpha', 0.5, ...       % overlay transparency
    'min_overlap', 0.2);           % require 20% of the track inside the region
```

`max_tracts` defaults to 2000 in whole-brain mode and to unlimited otherwise, so
setting it explicitly is usually worthwhile on dense track sets.

---

## Surveying every region

Grid mode renders all parcellation regions at once, and sequential mode steps
through them one at a time:

```matlab
visualizeTractography(run_dir, 'mode', 'grid');
visualizeTractography(run_dir, 'mode', 'grid', 'filter_mode', 'start_in');
visualizeTractography(run_dir, 'mode', 'sequential');
```

---

## Headless export

To export without an interactive session, use the shell wrapper, which renders
eight standard anatomical views through `visualizeTractographyAngles`:

```bash
# All tracks, eight views, into <run_dir>/figures/
./bin/run_visualization.sh hinec_runs/run_<TIMESTAMP>_hinec_dti/

# Restricted to regions 5, 10 and 15, at 600 dpi
./bin/run_visualization.sh hinec_runs/run_<TIMESTAMP>_hinec_dti/ png '5,10,15' 600
```

From MATLAB, `visualizeTractography` writes into a directory with
`export_dir` / `export_format` / `export_dpi`:

```matlab
visualizeTractography(run_dir, 'mode', 'region', 'region', 5, ...
    'export_dir', 'figures/', 'export_format', 'pdf', 'export_dpi', 600);
```

---

## Common use cases

| goal | call |
|---|---|
| everything related to a region | `'mode','region','region',id` |
| outgoing pathways | `..., 'filter_mode','start_in'` |
| incoming pathways | `..., 'filter_mode','end_in'` |
| inter-region connections | `..., 'filter_mode','connect_to'` |
| tract organisation | `..., 'color_mode','fa'` |

---

## Troubleshooting

**"Region X does not exist"** — the ID is not present in
`nim.parcellation_mask`. Check with `unique(nim.parcellation_mask(:))`, and
confirm the parcellation step ran.

**No tracks found for a region** — lower `min_overlap` or
`min_region_points`, switch from `start_in` to `pass_through`, or check whether
the region lies outside white matter. Small atlas labels are a common cause: the
JHU parcellation is thin (about 9,600 labelled voxels in the ISMRM data, with
individual labels down to a couple of dozen), so many regions intersect very few
streamlines. `nim_roi_mask` accepts a dilation argument for this reason.

**Tracks file not found** — pass the run directory and let auto-detection find
it, or check the exact name under `<run_dir>/tractography/`.

**Slow or unreadable rendering** — lower `max_tracts`, use
`color_mode: 'uniform'`, and set `show_region` false.
