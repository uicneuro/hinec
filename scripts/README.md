# Supported scripts

`bin/` contains the main pipeline launchers. The scripts here support those
launchers or provide documented commands for working with their outputs.

| Script | Use |
|---|---|
| `hinec_to_trk.py` | Convert saved MATLAB tracks to TRK; called by `bin/run_ismrm_scoring.sh`. |
| `mrtrix_peaks_to_hinec.py` | Convert MRtrix `sh2peaks` NIfTI output to a CSD peak MAT file on a HINEC reference grid. Use it with `tractography.csd.peaks_file`. |
| `ismrm_report_scores.py` | Print a whole-brain or seeded-bundle score; called by `bin/run_ismrm_scoring.sh`. |
| `build_ismrm_scoring_config.py` | Prepare the merged ISMRM scorer config. |
| `compare_ismrm_results.py` | Compare scores from completed runs. |
| `validate_ismrm_tractography.py` | Run standalone ISMRM validation diagnostics. |
| `FastTractographyViewer.py` | View a generated slice cache; launched by `bin/viewSlices.sh`. |
| `tractography_slice_gui.py` | Open the interactive MATLAB slice-viewer controls. |

Private experiments and one-off diagnostics are kept outside this directory.

To use a peak field reconstructed with MRtrix on the same DWI, convert its
`sh2peaks` NIfTI image to HINEC's voxel coordinates:

```bash
python scripts/mrtrix_peaks_to_hinec.py peaks.nii.gz dwi.nii.gz peaks.mat
```

Then set `tractography.field: csd` and `tractography.csd.peaks_file: peaks.mat`
in the run configuration. The MAT file must match the DWI grid; HINEC checks
its dimensions and active peak values before tracking.
