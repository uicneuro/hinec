# Distributed Workflow Guide

## Overview

The Fast Tractography Slice Viewer separates rendering from viewing, so that the
machine holding the data and MATLAB need not be the machine you browse the results
on. Three steps:

1. **Server** (high resources) — render every slice image with MATLAB
2. **Transfer** — copy the slice cache to the local computer
3. **Local** (no MATLAB) — browse the slices in the Python GUI

---

## The cache

The **cache** is a directory containing:

- **Pre-generated slice images** (PNG/JPG files)
- **Metadata** (JSON files with configuration)
- **NO original .mat files** (they stay on server)

```
cache/
├── datasets/
│   └── {dataset_hash}/
│       ├── metadata.json          ← Dataset information
│       └── parameters/
│           └── {param_hash}/
│               ├── config.json    ← Rendering parameters
│               ├── axial/
│               │   ├── axial_001.png
│               │   ├── axial_002.png
│               │   └── ... (all axial slices)
│               ├── sagittal/
│               │   └── ... (all sagittal slices)
│               └── coronal/
│                   └── ... (all coronal slices)
```

**Size**: Typically 50-500 MB depending on resolution and slice count.

---

## Step 1: Generate Slices on Server

### Requirements on Server
- MATLAB with Image Processing Toolbox
- HINEC pipeline installed
- Original data: `tracks.mat` + `nim.mat`
- Adequate storage (500 MB - 5 GB for cache)

### Run Generation

```matlab
% Put the HINEC sources on the path
addpath(genpath('src'));

% Basic usage (default settings)
cache_dir = generateSlices('tracks.mat', 'nim.mat');

% Custom output location
cache_dir = generateSlices('tracks.mat', 'nim.mat', '/export/slices');

% With custom options
options = struct();
options.image_resolution = [2048, 2048];  % Higher resolution
options.image_format = 'jpg';              % Smaller file size
options.tolerance = 3;                     % Thicker slices
cache_dir = generateSlices('tracks.mat', 'nim.mat', '/export/slices', options);
```

### Generation Options

| Option | Description | Default |
|--------|-------------|---------|
| `tolerance` | Slice thickness (voxels) | 2 |
| `image_resolution` | [width, height] pixels | [1024, 1024] |
| `image_format` | 'png' or 'jpg' | 'png' |
| `compression_level` | PNG: 0-9, JPG: 0-100 | 6 (PNG) / 90 (JPG) |
| `color_mode` | 'direction' or 'uniform' | 'direction' |
| `show_anatomy` | Show FA background | true |
| `alpha` | FA transparency | 0.6 |
| `parallel_workers` | Number of parallel workers | auto |

### Generation steps
1. Loads tractography data from .mat files
2. Generates slices for all three orientations (axial, sagittal, coronal)
3. Renders each slice as PNG/JPG image
4. Saves metadata and configuration
5. Validates completeness
6. Displays transfer instructions

**Time**: 5-30 minutes depending on data size and resolution.

---

## Step 2: Transfer Cache to Local Computer

### Option A: Compress and Transfer

```bash
# On server: Compress cache
cd /path/to/output
tar -czf slices.tar.gz cache/

# Transfer to local (choose one method)
scp slices.tar.gz user@local:/destination/
rsync -avz slices.tar.gz user@local:/destination/

# On local: Extract
tar -xzf slices.tar.gz
```

### Option B: Direct Sync

```bash
# From server to local
rsync -avz --progress server:/export/slices/ ~/local/slices/
```

### Option C: Physical Media

1. Copy cache directory to USB drive or external storage
2. Connect to local computer
3. Copy to local directory

### Verify Transfer

Check that the transferred directory has this structure:
```
slices/
└── datasets/
    └── {hash}/
        ├── metadata.json
        └── parameters/
```

---

## Step 3: View on Local Computer

### Requirements on Local Computer
- **No MATLAB required**
- Python 3.7+ with `tkinter`
- `scripts/FastTractographyViewer.py`
- Python packages: `pillow`, `numpy`

### Install Python Dependencies

```bash
pip install pillow numpy
```

### Launch Viewer

**Direct invocation**
```bash
python3 scripts/FastTractographyViewer.py /path/to/slices
```

**Launcher script** — `bin/viewSlices.sh` checks for Python and the `tkinter`,
`PIL` and `numpy` packages, installs any that are missing, validates the cache
directory, and then starts the viewer. Called with no argument it prompts for the
cache path.

```bash
./bin/viewSlices.sh /path/to/slices
```

### The viewer

The GUI opens with:

- Three orthogonal views (axial, sagittal, coronal), each with its own slider
- Sub-100 ms slice transitions, served from the pre-rendered images
- A "Synchronize crosshairs" toggle, in the control panel and under the View menu
- Performance monitoring in the status bar
- **File → Export Current Views…** to write the three visible slices out

### Keyboard shortcuts

Shortcuts act on the focused slice view.

| Key | Action |
|---|---|
| `←` / `→` | Move one slice |
| `↓` / `↑` | Move ten slices |

---

## Performance Characteristics

| Metric | Value |
|--------|-------|
| **Server Generation** | 5-30 minutes (one-time) |
| **Cache Size** | 50-500 MB typical |
| **Transfer Time** | Depends on network/media |
| **Local Navigation** | <100ms per slice |
| **Local Memory** | <500 MB |

---

## Workflow Comparison

### Traditional Approach
```
Local Computer: Open MATLAB → Load data → Render slice (5-30s delay!)
User changes slice → Wait 5-30s → Render next slice → Repeat...
```
**Problem**: Every slice change requires full computation.

### Fast Viewer Approach
```
Server: Generate all slices once (5-30 min) → Transfer to local
Local Computer: Open viewer → Instant navigation (<100ms)
User changes slice → Instant display → Smooth experience
```
**Benefit**: Pre-computation enables instant navigation.

---

## Troubleshooting

### Server Side

**Problem**: "Tracks file not found"
```matlab
% Check file paths
exist('tracks.mat', 'file')  % Should return 2
exist('nim.mat', 'file')     % Should return 2
```

**Problem**: "Out of memory"
```matlab
% Reduce resolution or use JPG
options.image_resolution = [512, 512];
options.image_format = 'jpg';
generateSlices('tracks.mat', 'nim.mat', '/export', options);
```

### Transfer

**Problem**: Transfer interrupted
```bash
# Resume with rsync
rsync -avz --partial --progress server:/export/slices/ ~/local/slices/
```

**Problem**: Permissions issue
```bash
# Fix permissions after transfer
chmod -R u+rwX ~/local/slices/
```

### Local Side

**Problem**: "Python not found"
```bash
# Check Python installation
python3 --version
which python3
```

**Problem**: "tkinter not found"
```bash
# On Ubuntu/Debian
sudo apt-get install python3-tk

# On macOS (usually built-in)
# Reinstall Python from python.org if needed
```

**Problem**: "Invalid cache directory"

- Make sure you're pointing to the ROOT cache directory
- Should contain `datasets/` subdirectory
- Check that transfer completed successfully

**Problem**: "Images not loading"

- Verify cache structure is intact
- Check file permissions: `ls -la /path/to/slices/datasets/`
- Re-extract if using compressed transfer

---

## Advanced Usage

### Multiple Datasets

Generate multiple datasets with different parameters:

```matlab
addpath(genpath('src'));

% High-resolution version
opts_hires = struct('image_resolution', [2048, 2048]);
generateSlices('tracks.mat', 'nim.mat', '/export/hires', opts_hires);

% Low-resolution version (faster transfer)
opts_lores = struct('image_resolution', [512, 512], 'image_format', 'jpg');
generateSlices('tracks.mat', 'nim.mat', '/export/lores', opts_lores);
```

### Batch Processing

Process multiple datasets:

```matlab
addpath(genpath('src'));

datasets = {'patient01', 'patient02', 'patient03'};

for i = 1:length(datasets)
    patient = datasets{i};
    tracks = sprintf('%s_tracks.mat', patient);
    nim = sprintf('%s_nim.mat', patient);
    output = sprintf('/export/%s', patient);

    fprintf('Processing %s...\n', patient);
    generateSlices(tracks, nim, output);
end
```

### Automated Transfer

Set up automated sync:

```bash
#!/bin/bash
# sync_slices.sh - Run on local computer

SERVER="user@server.edu"
REMOTE_DIR="/export/slices"
LOCAL_DIR="$HOME/tractography_data"

while true; do
    rsync -avz --progress "$SERVER:$REMOTE_DIR/" "$LOCAL_DIR/"
    echo "Sync complete. Waiting 1 hour..."
    sleep 3600
done
```

---

## Summary

✓ **Server**: `addpath(genpath('src')); generateSlices('tracks.mat', 'nim.mat', '/export/slices')`

✓ **Transfer**: `rsync -avz server:/export/slices/ ~/local/slices/`

✓ **Local**: `./bin/viewSlices.sh ~/local/slices/`

**Result**: Instant slice navigation with no MATLAB on local computer!