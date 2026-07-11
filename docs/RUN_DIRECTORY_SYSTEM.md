# HINEC Run Directory Organization System

## Overview

The HINEC pipeline now automatically organizes all generated files into timestamped run directories, preventing workspace clutter and enabling easy comparison between runs.

## Directory Structure

Each pipeline run creates an isolated directory:

```
hinec_runs/
├── run_20250118_143045_hinec_dti_cubic/     # Timestamped run directory
│   ├── config.yml                           # Copy of config used
│   ├── run_info.txt                         # Run metadata
│   ├── logs/
│   │   └── (MATLAB stdout/stderr)
│   ├── intermediate/
│   │   ├── sample.nii.gz                    # Preprocessed DWI
│   │   ├── sample_M.nii.gz                  # Brain mask
│   │   ├── sample_WM.nii.gz                 # White matter mask
│   │   ├── sample_GM.nii.gz                 # Gray matter mask
│   │   ├── sample_CSF.nii.gz                # CSF mask
│   │   └── parcellation_mask.nii.gz         # Parcellation
│   ├── output/
│   │   └── processed.mat                    # Final processed nim structure
│   └── tractography/
│       ├── tracks_hinec_2025-01-18_14_30_45.mat
│       └── diagnostics/
│           └── track_statistics.txt
│
├── run_20250118_150230_hinec_dti_fast/
│   └── (same structure)
│
└── latest -> run_20250118_150230_hinec_dti_fast/  # Symlink to most recent
```

## Benefits

✅ **Clean Workspace**: No scattered files in project root

✅ **Reproducibility**: Full run provenance with copied config

✅ **Easy Comparison**: Side-by-side run analysis

✅ **Simple Cleanup**: Delete old runs without hunting for files

✅ **Safe Experimentation**: Runs don't interfere with each other

✅ **Audit Trail**: Complete history of experiments

## Usage

### Basic Usage

```bash
# Run with default config - automatic run directory
./run_hinec.sh sample sample_output.mat

# Run with custom config
./run_hinec.sh sample sample_output.mat config/hinec_dti_cubic.yml
```

### Run Directory Contents

**config.yml**: Exact copy of configuration used for reproducibility

**run_info.txt**: Complete run metadata

- Run ID and timestamp
- Configuration used
- System information (MATLAB version, platform)
- Git commit hash (if available)
- Directory structure

**logs/**: All pipeline output (currently in startup log, future versions will redirect here)

**intermediate/**: All preprocessing artifacts

- Preprocessed DWI
- Brain mask (improved)
- Tissue masks (WM, GM, CSF)
- Parcellation mask

**output/**: Final processed data

- Main `.mat` file with complete nim structure
- Includes `nim.run_info` field for traceability

**tractography/**: Tractography results

- Track files (`.mat` format)
- **diagnostics/track_statistics.txt**: Detailed statistics
  - Track counts and lengths
  - Seeding statistics
  - Computation time
  - Parameters used

### Accessing Results

**Latest run** (via symlink):
```bash
ls -l hinec_runs/latest/
matlab -batch "load('hinec_runs/latest/output/processed.mat')"
```

**List all runs**:
```bash
ls -lt hinec_runs/
```

**Compare configs**:
```bash
diff hinec_runs/run_20250118_143045_hinec_dti_cubic/config.yml \
     hinec_runs/run_20250118_150230_hinec_dti_fast/config.yml
```

**Compare track counts**:
```bash
grep "Total tracks" hinec_runs/run_*/tractography/diagnostics/track_statistics.txt
```

**View run metadata**:
```bash
cat hinec_runs/run_20250118_143045_hinec_dti_cubic/run_info.txt
```

## File Naming Convention

### Run Directory Names

Format: `run_YYYYMMDD_HHMMSS_<config_preset>`

Examples:

- `run_20250118_143045_hinec_dti_cubic`
- `run_20250118_150230_hinec_dti_fast`
- `run_20250118_162015_irontract`
- `run_20250118_173200_hinec_default`

The config preset name is extracted from the YAML filename.

### Track Files

Format: `tracks_<algorithm>_YYYY-MM-DD_HH_MM_SS.mat`

Examples:

- `tracks_hinec_2025-01-18_14_30_45.mat`
- `tracks_standard_2025-01-18_15_02_30.mat`

## Cleanup and Management

### Delete Old Runs

```bash
# Delete specific run
rm -rf hinec_runs/run_20250118_143045_hinec_dti_cubic/

# Delete all runs before specific date
find hinec_runs/ -name "run_202501*" -type d -exec rm -rf {} +

# Keep only last 5 runs
ls -t hinec_runs/ | tail -n +6 | xargs -I {} rm -rf hinec_runs/{}
```

### Archive Runs

```bash
# Archive a run for long-term storage
tar -czf run_20250118_143045_hinec_dti_cubic.tar.gz \
    hinec_runs/run_20250118_143045_hinec_dti_cubic/

# Extract archived run
tar -xzf run_20250118_143045_hinec_dti_cubic.tar.gz -C hinec_runs/
```

## Integration with MATLAB

### Direct MATLAB Usage

You can also use the run directory system directly from MATLAB:

```matlab
% Load configuration
config = load_config_yaml('config/hinec_dti_cubic.yml');

% Create run directory
run_info = create_run_directory('config/hinec_dti_cubic.yml');

% Run preprocessing with run directory
main('sample', 'processed.mat', config, run_info);

% Run tractography with run directory
output_mat = fullfile(run_info.output_dir, 'processed.mat');
runTractography(output_mat, config, run_info);

% Access run directory
fprintf('Results in: %s\n', run_info.run_dir);
```

### Custom Run Names

```matlab
% Create run with custom description
run_info = create_run_directory('config/irontract.yml', ...
    'run_name', 'my_experiment_v2', ...
    'description', 'Testing new seeding strategy');
```

## Troubleshooting

### Issue: Can't find latest results

**Solution**: Use the `latest` symlink
```bash
ls hinec_runs/latest/
```

### Issue: Run directory not created

**Check**: MATLAB path includes `nim_utils/`
```matlab
addpath('nim_utils');
```

### Issue: Want to disable run directory system

**Solution**: Currently automatic. Legacy mode (files in workspace) will be added in future version if needed.

### Issue: Symlink not working on Windows

**Note**: On Windows, the system creates a text file instead of a symlink. Read the file to get the path:
```bash
cat hinec_runs/latest  # Shows path to latest run directory
```

## Best Practices

1. **Keep Configs**: Never delete the YAML config files - they're copied to run directories for reproducibility

2. **Document Experiments**: Use run_info.txt or create additional README files in run directories

3. **Regular Cleanup**: Periodically archive or delete old experimental runs to save disk space

4. **Comparison Workflow**:
   ```bash
   # Quick comparison of two runs
   diff hinec_runs/run_A/config.yml hinec_runs/run_B/config.yml
   cat hinec_runs/run_A/tractography/diagnostics/track_statistics.txt
   cat hinec_runs/run_B/tractography/diagnostics/track_statistics.txt
   ```

5. **Version Control**: The run system works well with git - run directories are automatically excluded (add to .gitignore)

## Technical Details

### Run Directory Creation

The `create_run_directory.m` function:

1. Creates timestamped directory name from config file
2. Creates all subdirectories
3. Copies configuration file
4. Generates run metadata with system info and git commit
5. Updates `latest` symlink
6. Returns `run_info` structure with all paths

### Integration Points

**main.m**:

- Accepts optional `run_info` parameter
- Redirects all intermediate files to `run_info.intermediate_dir`
- Saves output to `run_info.output_dir`
- Stores `run_info` in nim structure

**runTractography.m**:

- Accepts optional `run_info` parameter
- Saves tracks to `run_info.tractography_dir`
- Generates diagnostics in `run_info.diagnostics_dir`

**run_hinec.sh**:

- Calls `create_run_directory()` before pipeline
- Passes `run_info` to both main() and runTractography()
- Updates progress messages with run directory location

### Backward Compatibility

The system is fully backward compatible:

- Old scripts without run_info parameter still work
- Files fall back to current directory behavior
- YAML system works independently of run directories

## Future Enhancements

Planned features:

- Optional flag to disable run directory system
- Automatic cleanup of runs older than N days
- Run comparison tool
- Integration with visualization pipeline
