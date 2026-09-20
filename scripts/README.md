# Supported scripts

`bin/` contains the main pipeline launchers. The scripts here support those
launchers or provide documented commands for working with their outputs.

| Script | Use |
|---|---|
| `hinec_to_trk.py` | Convert saved MATLAB tracks to TRK; called by `bin/run_ismrm_scoring.sh`. |
| `ismrm_report_scores.py` | Print a whole-brain or seeded-bundle score; called by `bin/run_ismrm_scoring.sh`. |
| `build_ismrm_scoring_config.py` | Prepare the merged ISMRM scorer config. |
| `compare_ismrm_results.py` | Compare scores from completed runs. |
| `validate_ismrm_tractography.py` | Run standalone ISMRM validation diagnostics. |
| `FastTractographyViewer.py` | View a generated slice cache; launched by `bin/viewSlices.sh`. |
| `tractography_slice_gui.py` | Open the interactive MATLAB slice-viewer controls. |
| `generate_presentation_figures.m` | Generate figures using `src/nim_presentation/`. |

Private experiments and one-off diagnostics are kept outside this directory.
