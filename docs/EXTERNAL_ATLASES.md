# Using an external atlas

The full `bin/run_hinec.sh` pipeline accepts a user-supplied label atlas through
`preprocessing.atlas_file`. This replaces built-in atlas registration in both
preprocessing and the final parcellation stage. It works with human or nonhuman
anatomy, provided the atlas has already been aligned to the DWI image grid.

```yaml
preprocessing:
  mask_file: data/my_subject/brain_mask_in_dwi.nii.gz
  atlas_file: data/my_subject/atlas_in_dwi.nii.gz
  atlas_labels_file: data/my_subject/atlas_labels.tsv
  use_t1_registration: false
  register_to_mni: false

tractography:
  seeding:
    roi: [Callosal subset]
```

```bash
./bin/run_hinec.sh data/my_subject/subject subject.mat config/my_subject.yml
```

Paths are resolved from the pipeline working directory, just like `mask_file`;
absolute paths are also accepted. CLI overrides use the same keys. The launcher
and raw-DWI tensor estimation workflow are unchanged.

## Brain masks and atlases have different roles

- `mask_file` is the brain/processing mask; it does not define named regions.
- `atlas_file` is a 3-D integer parcellation for region labels, seeding, and
  anatomical selection. It does not replace the brain mask.
- `bundle_roi_dir` can subsequently replace the active parcellation with bundle
  masks, preserving the external atlas under `parcellation_mask_external` and
  its names under `atlas_labels_external`.

When `atlas_file` is empty, built-in `atlas_type` behavior remains available.
When supplied, `atlas_file` takes precedence over `atlas_type`; there is no
fallback to a human atlas if validation fails. T1 availability does not override
an explicit `use_t1_registration: false`. With an external atlas, T1-to-MNI
registration runs only when `register_to_mni: true` is requested. For nonhuman
data, keep human MNI registration disabled.

## Grid and label requirements

The atlas must be `.nii` or `.nii.gz`, with the same spatial dimensions, units,
voxel spacing, and affine as the DWI. The pipeline does not register, transpose,
or resize supplied atlases. If your atlas is in template or T1 space, first
register it to the subject's DWI space using the appropriate species-specific
reference and nearest-neighbour interpolation. Matching headers are necessary
but do not prove anatomical alignment; inspect an overlay before tracking.

Voxel values must be finite nonnegative integers within int32 range. Zero is
background. Probability maps, negative values, and all-background volumes are
rejected. Label values are preserved exactly, without the index offsets used
by some atlas XML formats. Atlas alignment is checked against the raw DWI before
expensive preprocessing and against the processed DWI before final attachment.

`atlas_labels_file` is optional. It is a UTF-8 tab-separated text file with these
columns (the gaps below are tabs):

```text
index	name
0	Background
7	Callosal subset
42	Other region
```

Every positive label present in the image must appear in a supplied table.
Unused table entries and an optional background row are allowed. IDs must be
unique nonnegative integers, and names must be nonempty and unique ignoring
case. Without a table, names are generated as `Region_7`, `Region_42`, etc.;
numeric ROI selections such as `roi: [7]` also work. A labels table without an
`atlas_file` is rejected.

## Saved outputs and reuse

The run's `intermediate/parcellation_mask.nii.gz` contains the staged atlas.
`external_atlas_labels.mat` records its names, source paths, and coordinate
policy. The processed `nim` also embeds the map and sets `atlas_type` to
`external`, so ROI selection does not depend on discovering a human XML file.
Existing XML sidecars cannot overwrite explicit external or bundle names.

Use a new output MAT path when applying an external atlas. The full pipeline
rejects an existing output MAT instead of silently reusing its old anatomy.
Existing preprocessed DWI may still be reused; the supplied atlas is validated
and attached again. To rerun raw preprocessing, use fresh writable input staging
without the cached `<prefix>.nii.gz` output.
