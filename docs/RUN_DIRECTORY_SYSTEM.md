# HINEC Run Directory Organization System

## Overview

Every pipeline invocation writes into its own timestamped directory under `hinec_runs/`,
carrying the config it ran with. Two properties follow from that: nothing lands in the project
root, and a finished run states its own parameters, so results stay attributable after the
fact.

The division of labour is worth stating up front, because it is easy to get backwards. The
processed `nim`, the preprocessed references and the CSD cache belong to the **data layer**
(next to the input, e.g. `data/ismrm2015/`). Run directories hold **tractography outputs**.
`main.m` writes the canonical `nim` to the data layer and a copy into the run directory for
that run's own provenance; never treat a run directory as the home of a `nim`.

## Directory Structure

```
hinec_runs/
├── run_<timestamp>_<config>/                # e.g. run_20250118_143045_hinec_dti_cubic
│   ├── config.yml                           # Copy of the config used
│   ├── run_info.txt                         # Run metadata
│   ├── overrides.txt                        # --set overrides, when any were given
│   ├── logs/
│   │   └── pipeline.log                     # Live MATLAB stdout/stderr
│   ├── intermediate/
│   │   ├── sample.nii.gz                    # Preprocessed DWI
│   │   ├── sample_M.nii.gz                  # Brain mask
│   │   ├── sample_WM_mask.nii.gz            # White matter mask
│   │   ├── sample_GM_mask.nii.gz            # Gray matter mask
│   │   ├── sample_CSF_mask.nii.gz           # CSF mask
│   │   ├── parcellation_mask.nii.gz         # Parcellation
│   │   └── SOURCE.txt                       # (run_tractography.sh) frozen DWI reference
│   ├── output/
│   │   └── processed.mat                    # Copy of the processed nim, with run_info
│   ├── tractography/
│   │   ├── tracks_hinec_<timestamp>.mat
│   │   └── diagnostics/
│   │       └── track_statistics.txt
│   ├── scoring/                             # (--score) Renauld-2023 scorer output
│   └── figures/                             # (run_visualization.sh) exported figures
│
└── latest -> run_<timestamp>_<config>/      # Symlink to the most recent run
```

`run_tractography.sh` copies — rather than symlinks — the DWI reference into its own
`intermediate/` and records where it came from in `SOURCE.txt`. Freezing the reference means
a later reprocess of the source data cannot silently change what an old run was scored
against.

## Usage

### Basic Usage

```bash
# Full pipeline; the run directory is created for you
./bin/run_hinec.sh data/parcellation_sample/sample sample_output.mat

# With a specific config
./bin/run_hinec.sh data/parcellation_sample/sample sample_output.mat config/hinec_dti_cubic.yml

# Tractography only, against an already-preprocessed nim
./bin/run_tractography.sh hinec_dti_cubic --score
```

### Run Directory Contents

**config.yml**: Exact copy of configuration used for reproducibility

**run_info.txt**: Complete run metadata

- Run ID and timestamp
- Configuration used
- System information (MATLAB version, platform)
- Git commit hash (if available)
- Directory structure

**logs/pipeline.log**: the live MATLAB stdout/stderr for the run; tail it to watch progress

**overrides.txt**: the `--set key=value` overrides applied on the command line, one per line.
Written only when overrides were given, so its presence tells you the run's `config.yml` is
not the whole story.

**intermediate/**: preprocessing artifacts

- Preprocessed DWI
- Brain mask (improved)
- Tissue masks (WM, GM, CSF)
- Parcellation mask

**output/**: a copy of the processed `nim`, with `nim.run_info` embedded for traceability. The
canonical copy lives in the data directory.

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
diff hinec_runs/run_<A>/config.yml hinec_runs/run_<B>/config.yml
cat hinec_runs/run_*/overrides.txt
```

**Compare track counts**:
```bash
grep "Total tracks" hinec_runs/run_*/tractography/diagnostics/track_statistics.txt
```

**View run metadata**:
```bash
cat hinec_runs/latest/run_info.txt
```

## File Naming Convention

### Run Directory Names

Format: `run_YYYYMMDD_HHMMSS_<config_preset>`

For example, a run launched with `config/hinec_dti_cubic.yml` at 14:30:45 on 18 January 2025
becomes `run_20250118_143045_hinec_dti_cubic`. The preset name is the config file's basename.

`run_tractography.sh` additionally tags the name with any `--set` overrides, so a sweep
produces distinguishable directories rather than a row of identical timestamps.

### Track Files

Format: `tracks_<algorithm>_YYYY-MM-DD_HH_MM_SS.mat`

One file per run, named for the tracker that produced it: `tracks_standard_*`,
`tracks_hinec_*` or `tracks_mmf_*`. For example,
`tracks_hinec_2025-01-18_14_30_45.mat`.

## Cleanup and Management

### Delete Old Runs

```bash
# Delete one run
rm -rf hinec_runs/run_<timestamp>_<config>/

# Delete every run from a given month
find hinec_runs/ -maxdepth 1 -name "run_202501*" -type d -exec rm -rf {} +

# Keep only the 5 most recent
ls -t hinec_runs/ | tail -n +6 | xargs -I {} rm -rf hinec_runs/{}
```

### Archive Runs

```bash
# Archive a run for long-term storage
run=run_<timestamp>_<config>
tar -czf "${run}.tar.gz" "hinec_runs/${run}/"

# Extract it again
tar -xzf "${run}.tar.gz" -C hinec_runs/
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

**Check**: the MATLAB path includes `src/nim_utils/`, where `create_run_directory` lives.
```matlab
addpath(genpath('.'));   % what the shell launchers do
```

### Issue: Running without a run directory

Pass no `run_info` and both `main` and `runTractography` fall back to writing beside their
inputs — `main` next to the image path, `runTractography` into `tractography_results/`. Every
code path that touches outputs branches on whether a `run_info` was supplied.

### Issue: `latest` is a file, not a symlink

On Windows, `create_run_directory` writes a text file containing the path instead of creating
a symlink. Read the file to get the path:
```bash
cat hinec_runs/latest
```
A failure to create the link is a warning, not an error; the run directory itself is still
complete.

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

5. **Version Control**: `hinec_runs/` is already in `.gitignore`, so runs never enter git.

## Technical Details

### Run Directory Creation

`create_run_directory.m`:

1. Builds a timestamped directory name from the config file's basename
2. Creates the subdirectories (`mkdir` is idempotent, so reusing a name is safe)
3. Copies the config file to `config.yml`
4. Writes `run_info.txt` with the run ID, timestamp, MATLAB version, platform, working
   directory, and the git commit, branch and clean/modified status when git is available
5. Updates the `latest` link
6. Returns the `run_info` struct

The returned struct carries the paths every downstream stage writes through:

| Field | Contents |
|---|---|
| `.run_dir` | The run directory itself |
| `.config_file` | Path to the copied `config.yml` |
| `.logs_dir` | `<run_dir>/logs` |
| `.intermediate_dir` | `<run_dir>/intermediate` |
| `.output_dir` | `<run_dir>/output` |
| `.tractography_dir` | `<run_dir>/tractography` |
| `.diagnostics_dir` | `<run_dir>/tractography/diagnostics` |
| `.run_id` | The directory's own name |
| `.timestamp` | `yyyymmdd_HHMMSS` |

Optional name-value arguments: `base_dir` (default `hinec_runs`), `run_name` (overrides the
generated name entirely) and `description` (recorded in `run_info.txt`).

### Integration Points

**main.m** takes `run_info` as its fourth positional argument, redirects intermediate files to
`run_info.intermediate_dir`, writes its output copy to `run_info.output_dir`, and stores
`run_info` on the `nim`.

**runTractography.m** takes `run_info` as its third positional argument, writes tracks to
`run_info.tractography_dir` and diagnostics to `run_info.diagnostics_dir`.

**bin/run_hinec.sh** creates the run directory before anything else, then passes the same
`run_info` to both `main` and `runTractography` in one `matlab -batch` process.

**bin/run_tractography.sh** creates its own run directory, freezes the DWI reference into it,
and passes `run_info` to `runTractography` only.

### Backward Compatibility

Callers that pass no `run_info` still work: `main` and `runTractography` fall back to writing
beside their inputs. The YAML config system is independent of run directories and works either
way.
