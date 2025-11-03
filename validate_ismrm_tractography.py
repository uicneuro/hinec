#!/usr/bin/env python3
"""
validate_ismrm_tractography.py - ISMRM Tractography Validation & Scoring

Comprehensive validation of HINEC tractography against ISMRM 2015 ground truth.
Converts MATLAB tracks, performs ROI-based segmentation, and generates accuracy report.

Usage:
    python validate_ismrm_tractography.py \\
        --nim-file data.mat \\
        --tracks-file tractography_results/tracks_standard.mat \\
        --scoring-dir ISMRM_data/scoring_data_Renauld2023 \\
        --output-dir validation_results/

Metrics Computed:
    - Bundle coverage: % of ground truth covered by submission
    - Bundle overreach: % of submission not in ground truth
    - Valid bundles: % of bundles passing anatomical checks
    - Invalid connections: Anatomically impossible streamlines
    - Spatial overlap: Weighted Dice coefficient
    - Track statistics: Count, length, coverage per bundle

Output:
    - HTML report with visualizations
    - Per-bundle metrics (CSV)
    - Segmented bundles (TRK files)
    - Coordinate alignment diagnostics
"""

import argparse
import json
import os
import sys
from pathlib import Path
from datetime import datetime
import numpy as np

# Check dependencies
try:
    import nibabel as nib
    from nibabel.streamlines import Tractogram, load as load_trk, save as save_trk
except ImportError:
    print("Error: nibabel not installed. Install with: pip install nibabel")
    sys.exit(1)

try:
    import scipy.io
    from scipy.spatial.distance import cdist
except ImportError:
    print("Error: scipy not installed. Install with: pip install scipy")
    sys.exit(1)

try:
    import h5py
    H5PY_AVAILABLE = True
except ImportError:
    print("Warning: h5py not installed. MATLAB v7.3 files will not be readable.")
    print("Install with: pip install h5py")
    H5PY_AVAILABLE = False

try:
    from dipy.tracking.streamline import set_number_of_points
    from dipy.tracking.distances import bundles_distances_mam
except ImportError:
    print("Warning: dipy not installed. Some metrics will be unavailable.")
    print("Install with: pip install dipy")
    DIPY_AVAILABLE = False
else:
    DIPY_AVAILABLE = True


class ISMRMValidator:
    """ISMRM Tractography Validation Engine"""

    def __init__(self, scoring_dir, output_dir):
        self.scoring_dir = Path(scoring_dir)
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)

        # Load configurations
        self.seg_config = self._load_json(self.scoring_dir / 'config_file_segmentation.json')
        self.tract_config = self._load_json(self.scoring_dir / 'config_file_tractometry.json')

        # Results storage
        self.results = {
            'bundles': {},
            'summary': {},
            'diagnostics': {}
        }

    def _load_json(self, filepath):
        """Load JSON configuration file"""
        with open(filepath, 'r') as f:
            return json.load(f)

    def load_hinec_tracks(self, mat_file):
        """
        Load HINEC tractography from MATLAB .mat file

        Supports both standard and v7.3 (HDF5) formats.

        Returns:
            list of Nx3 numpy arrays (0-based voxel coordinates)
        """
        print(f"\n=== Loading HINEC Tractography ===")
        print(f"File: {mat_file}")

        # Try standard MATLAB format first
        try:
            mat_data = scipy.io.loadmat(mat_file)
            print("Format: MATLAB standard (.mat)")

            # Try different field names
            if 'tracks' in mat_data:
                tracks_cell = mat_data['tracks']
            elif 'tracts' in mat_data:
                tracks_cell = mat_data['tracts']
            else:
                raise ValueError("No 'tracks' or 'tracts' field in .mat file")

            # Convert MATLAB cell array to Python list
            tracks = []
            for i in range(tracks_cell.shape[0]):
                track = tracks_cell[i, 0]
                if track.size > 0:
                    # Convert 1-based MATLAB to 0-based Python
                    track_0based = track.astype(np.float32) - 1.0
                    tracks.append(track_0based)

        except NotImplementedError as e:
            # This is a v7.3 file, need h5py
            if "v7.3" in str(e) or "HDF" in str(e):
                if not H5PY_AVAILABLE:
                    raise ImportError(
                        "This is a MATLAB v7.3 file (HDF5 format).\n"
                        "Please install h5py: pip install h5py"
                    )

                print("Format: MATLAB v7.3 (HDF5)")
                tracks = self._load_hdf5_tracks(mat_file)
            else:
                raise

        print(f"✓ Loaded {len(tracks)} streamlines")

        # Coordinate statistics
        if len(tracks) > 0:
            all_coords = np.vstack(tracks)
            print(f"Coordinate range (0-based voxels):")
            print(f"  X: [{all_coords[:,0].min():.1f}, {all_coords[:,0].max():.1f}]")
            print(f"  Y: [{all_coords[:,1].min():.1f}, {all_coords[:,1].max():.1f}]")
            print(f"  Z: [{all_coords[:,2].min():.1f}, {all_coords[:,2].max():.1f}]")

            self.results['diagnostics']['input_coords'] = {
                'x_range': [float(all_coords[:,0].min()), float(all_coords[:,0].max())],
                'y_range': [float(all_coords[:,1].min()), float(all_coords[:,1].max())],
                'z_range': [float(all_coords[:,2].min()), float(all_coords[:,2].max())]
            }

        return tracks

    def _load_hdf5_tracks(self, mat_file):
        """
        Load tracks from MATLAB v7.3 (HDF5) file

        Args:
            mat_file: Path to .mat file

        Returns:
            list of Nx3 numpy arrays
        """
        with h5py.File(mat_file, 'r') as f:
            # Try different field names
            if 'tracks' in f:
                tracks_ref = f['tracks']
            elif 'tracts' in f:
                tracks_ref = f['tracts']
            else:
                raise ValueError("No 'tracks' or 'tracts' field in HDF5 file")

            print(f"  HDF5 structure: {list(f.keys())}")

            # v7.3 stores cell arrays as object references
            tracks = []

            # Check if it's a cell array (references) or direct array
            if tracks_ref.dtype == np.dtype('O'):
                # Cell array format - each element is a reference
                # Shape can be (N, 1) column vector or (1, N) row vector
                num_tracks = tracks_ref.shape[0] * tracks_ref.shape[1]
                print(f"  Loading {num_tracks} track references...")

                # Flatten to handle both (N,1) and (1,N) shapes
                # Need to read into numpy array first to use flatten()
                tracks_array = tracks_ref[...]  # Read entire dataset
                tracks_flat = tracks_array.flatten()

                for i in range(num_tracks):
                    # Dereference each track
                    track_ref = tracks_flat[i]
                    track = f[track_ref][()]  # Load the actual data

                    if track.size > 0:
                        # HDF5 stores in column-major (Fortran) order
                        # MATLAB stores as 3xN, we need Nx3
                        if track.shape[0] == 3:
                            track = track.T  # Transpose to Nx3
                        elif track.shape[1] == 3:
                            pass  # Already Nx3
                        else:
                            print(f"    Warning: Unexpected track shape {track.shape}, skipping")
                            continue

                        # Convert 1-based MATLAB to 0-based Python
                        track_0based = track.astype(np.float32) - 1.0
                        tracks.append(track_0based)

                    if (i + 1) % 10000 == 0:
                        print(f"    Loaded {i + 1} / {num_tracks} tracks...")

            else:
                # Direct array format (less common)
                print("  Loading direct array format...")
                track_data = tracks_ref[()]

                # Split into individual tracks if needed
                # This depends on how MATLAB saved it - may need adjustment
                if track_data.ndim == 2 and track_data.shape[1] == 3:
                    # Single track or concatenated tracks
                    track_0based = track_data.astype(np.float32) - 1.0
                    tracks.append(track_0based)
                else:
                    raise ValueError(f"Unexpected track data shape: {track_data.shape}")

        return tracks

    def convert_to_world_space(self, tracks, reference_nii):
        """
        Convert tracks from voxel to world coordinates

        Args:
            tracks: List of Nx3 arrays (0-based voxel coords)
            reference_nii: Path to reference NIfTI

        Returns:
            List of Nx3 arrays (world coords)
        """
        print(f"\n=== Converting to World Space ===")
        ref_img = nib.load(reference_nii)
        affine = ref_img.affine

        print(f"Reference: {reference_nii}")
        print(f"Dimensions: {ref_img.shape}")
        print(f"Voxel size: {ref_img.header.get_zooms()[:3]} mm")

        # Convert to world coordinates
        world_tracks = []
        for track in tracks:
            # Add homogeneous coordinate
            track_homog = np.column_stack([track, np.ones(len(track))])
            # Apply affine
            track_world = track_homog @ affine.T
            world_tracks.append(track_world[:, :3].astype(np.float32))

        # Statistics
        all_coords = np.vstack(world_tracks)
        print(f"World coordinate range:")
        print(f"  X: [{all_coords[:,0].min():.1f}, {all_coords[:,0].max():.1f}]")
        print(f"  Y: [{all_coords[:,1].min():.1f}, {all_coords[:,1].max():.1f}]")
        print(f"  Z: [{all_coords[:,2].min():.1f}, {all_coords[:,2].max():.1f}]")

        self.results['diagnostics']['world_coords'] = {
            'x_range': [float(all_coords[:,0].min()), float(all_coords[:,0].max())],
            'y_range': [float(all_coords[:,1].min()), float(all_coords[:,1].max())],
            'z_range': [float(all_coords[:,2].min()), float(all_coords[:,2].max())]
        }

        return world_tracks

    def load_ground_truth(self, bundle_name):
        """Load ground truth streamlines for a bundle"""
        # Check if bundle exists in tractometry config
        if bundle_name not in self.tract_config:
            # Bundle is in segmentation config but not tractometry (e.g., CC_temporal)
            # These are sub-bundles that don't have separate ground truth
            return None

        gt_path = self.scoring_dir / self.tract_config[bundle_name]['gt_mask']
        if not gt_path.exists():
            return None

        gt_tractogram = load_trk(str(gt_path), 'same')
        return list(gt_tractogram.streamlines)

    def segment_bundle_roi(self, tracks, bundle_name, bundle_config, dwi_to_t1_scale=None):
        """
        Segment bundle using ROI-based filtering

        Applies anatomical constraints from config:
        - all_mask: All points must pass through
        - any_mask: At least one point must pass through
        - head/tail: Endpoints must be in these regions
        - length constraints: Geometric filtering

        Args:
            tracks: List of Nx3 arrays (voxel coords in DWI space)
            bundle_name: Name of bundle
            bundle_config: ROI configuration
            dwi_to_t1_scale: Scaling factor from DWI to T1 space (default: 2.0 for 2mm→1mm)
        """
        print(f"\n  Segmenting {bundle_name}...")

        # Load masks (in T1 space)
        all_mask = None
        if 'all_mask' in bundle_config:
            all_mask_path = self.scoring_dir / bundle_config['all_mask']
            if all_mask_path.exists():
                all_mask = nib.load(str(all_mask_path)).get_fdata() > 0

        any_mask = None
        if 'any_mask' in bundle_config:
            any_mask_path = self.scoring_dir / bundle_config['any_mask']
            if any_mask_path.exists():
                any_mask = nib.load(str(any_mask_path)).get_fdata() > 0

        head_mask = None
        if 'head' in bundle_config:
            head_path = self.scoring_dir / bundle_config['head']
            if head_path.exists():
                head_mask = nib.load(str(head_path)).get_fdata() > 0

        tail_mask = None
        if 'tail' in bundle_config:
            tail_path = self.scoring_dir / bundle_config['tail']
            if tail_path.exists():
                tail_mask = nib.load(str(tail_path)).get_fdata() > 0

        # Determine scaling factor
        if dwi_to_t1_scale is None:
            # Auto-detect: ISMRM uses 2mm DWI → 1mm T1 (scale by 2)
            dwi_to_t1_scale = 2.0

        # Filter tracks
        filtered = []
        for track in tracks:
            # Transform from DWI space to T1 space
            # DWI @ 90×108×90 (2mm) → T1 @ 180×216×180 (1mm)
            track_t1 = track * dwi_to_t1_scale

            # Convert to voxel indices
            track_vox = np.round(track_t1).astype(int)

            # Check bounds
            if not self._in_bounds(track_vox, all_mask.shape if all_mask is not None else (100, 100, 100)):
                continue

            # All mask check
            if all_mask is not None:
                if not self._passes_all_mask(track_vox, all_mask):
                    continue

            # Any mask check
            if any_mask is not None:
                if not self._passes_any_mask(track_vox, any_mask):
                    continue

            # Endpoint checks
            if head_mask is not None and tail_mask is not None:
                if not self._check_endpoints(track_vox, head_mask, tail_mask):
                    continue

            # Length constraints
            if not self._check_length_constraints(track, bundle_config):
                continue

            filtered.append(track)

        print(f"    {len(filtered)} / {len(tracks)} streamlines passed ({100*len(filtered)/len(tracks):.1f}%)")
        return filtered

    def _in_bounds(self, track_vox, shape):
        """Check if track is within volume bounds"""
        return np.all((track_vox >= 0) & (track_vox < np.array(shape)))

    def _passes_all_mask(self, track_vox, mask):
        """Check if all points pass through mask"""
        for point in track_vox:
            if np.all(point >= 0) and np.all(point < np.array(mask.shape)):
                if mask[tuple(point)] == 0:
                    return False
        return True

    def _passes_any_mask(self, track_vox, mask):
        """Check if any point passes through mask"""
        for point in track_vox:
            if np.all(point >= 0) and np.all(point < np.array(mask.shape)):
                if mask[tuple(point)] > 0:
                    return True
        return False

    def _check_endpoints(self, track_vox, head_mask, tail_mask):
        """Check if endpoints are in correct regions"""
        start = track_vox[0]
        end = track_vox[-1]

        # Check if start in head OR tail (bidirectional)
        start_in_head = head_mask[tuple(start)] > 0 if self._in_bounds(start.reshape(1,-1), head_mask.shape) else False
        start_in_tail = tail_mask[tuple(start)] > 0 if self._in_bounds(start.reshape(1,-1), tail_mask.shape) else False

        end_in_head = head_mask[tuple(end)] > 0 if self._in_bounds(end.reshape(1,-1), head_mask.shape) else False
        end_in_tail = tail_mask[tuple(end)] > 0 if self._in_bounds(end.reshape(1,-1), tail_mask.shape) else False

        return (start_in_head and end_in_tail) or (start_in_tail and end_in_head)

    def _check_length_constraints(self, track, config):
        """Check geometric length constraints"""
        # Total length
        if 'length' in config:
            length = self._compute_length(track)
            if not (config['length'][0] <= length <= config['length'][1]):
                return False

        # Axis-specific constraints
        extent = np.ptp(track, axis=0)  # max - min per axis

        if 'length_x' in config:
            if not (config['length_x'][0] <= extent[0] <= config['length_x'][1]):
                return False

        if 'length_y' in config:
            if not (config['length_y'][0] <= extent[1] <= config['length_y'][1]):
                return False

        if 'length_z' in config:
            if not (config['length_z'][0] <= extent[2] <= config['length_z'][1]):
                return False

        if 'length_x_abs' in config:
            if not (config['length_x_abs'][0] <= abs(extent[0]) <= config['length_x_abs'][1]):
                return False

        return True

    def _compute_length(self, track):
        """Compute track length in mm"""
        if len(track) < 2:
            return 0
        diffs = np.diff(track, axis=0)
        return np.sum(np.sqrt(np.sum(diffs**2, axis=1)))

    def compute_bundle_metrics(self, user_bundle, gt_bundle, bundle_name):
        """
        Compute accuracy metrics for a bundle

        Metrics:
            - Coverage: % of ground truth covered
            - Overreach: % of user bundle not in ground truth
            - Dice coefficient: Spatial overlap
            - Track statistics: Count, length, density
        """
        metrics = {
            'name': bundle_name,
            'user_count': len(user_bundle),
            'gt_count': len(gt_bundle) if gt_bundle else 0,
            'coverage': 0.0,
            'overreach': 0.0,
            'dice': 0.0,
            'avg_length_user': 0.0,
            'avg_length_gt': 0.0,
            'valid': len(user_bundle) > 0
        }

        if len(user_bundle) == 0:
            return metrics

        # Compute average lengths (convert numpy types to Python float)
        metrics['avg_length_user'] = float(np.mean([self._compute_length(t) for t in user_bundle]))

        if gt_bundle and len(gt_bundle) > 0:
            metrics['avg_length_gt'] = float(np.mean([self._compute_length(t) for t in gt_bundle]))

            # Compute coverage and overreach using MDF (if Dipy available)
            if DIPY_AVAILABLE:
                try:
                    # Resample to same number of points
                    user_resampled = [set_number_of_points(t, 20) for t in user_bundle]
                    gt_resampled = [set_number_of_points(t, 20) for t in gt_bundle]

                    # Compute distances
                    distances = bundles_distances_mam(user_resampled, gt_resampled)

                    # Coverage: % of GT covered by user tracks (distance < threshold)
                    threshold = 10.0  # mm
                    covered_gt = np.sum(np.min(distances, axis=0) < threshold)
                    metrics['coverage'] = float(100 * covered_gt / len(gt_bundle))

                    # Overreach: % of user tracks not covering GT
                    not_covering = np.sum(np.min(distances, axis=1) >= threshold)
                    metrics['overreach'] = float(100 * not_covering / len(user_bundle))

                    # Dice-like metric (convert to Python float)
                    metrics['dice'] = float(2 * covered_gt / (len(user_bundle) + len(gt_bundle)))

                except Exception as e:
                    print(f"      Warning: Could not compute distances: {e}")

        return metrics

    def generate_report(self, nim_file, tracks_file):
        """Generate HTML validation report"""
        print(f"\n=== Generating Report ===")

        report_path = self.output_dir / 'validation_report.html'

        # Calculate summary statistics
        total_bundles = len(self.results['bundles'])
        valid_bundles = sum(1 for b in self.results['bundles'].values() if b['valid'])
        total_user_tracks = sum(b['user_count'] for b in self.results['bundles'].values())
        total_gt_tracks = sum(b['gt_count'] for b in self.results['bundles'].values() if b['gt_count'] > 0)

        avg_coverage = np.mean([b['coverage'] for b in self.results['bundles'].values() if b['gt_count'] > 0])
        avg_overreach = np.mean([b['overreach'] for b in self.results['bundles'].values() if b['user_count'] > 0])

        self.results['summary'] = {
            'total_bundles': total_bundles,
            'valid_bundles': valid_bundles,
            'total_user_tracks': total_user_tracks,
            'total_gt_tracks': total_gt_tracks,
            'avg_coverage': float(avg_coverage),
            'avg_overreach': float(avg_overreach),
            'timestamp': datetime.now().isoformat()
        }

        # Generate HTML
        html = f"""<!DOCTYPE html>
<html>
<head>
    <title>ISMRM Tractography Validation Report</title>
    <style>
        body {{ font-family: Arial, sans-serif; margin: 20px; background: #f5f5f5; }}
        .container {{ max-width: 1200px; margin: 0 auto; background: white; padding: 30px; box-shadow: 0 0 10px rgba(0,0,0,0.1); }}
        h1 {{ color: #333; border-bottom: 3px solid #4CAF50; padding-bottom: 10px; }}
        h2 {{ color: #555; margin-top: 30px; }}
        .summary {{ background: #e8f5e9; padding: 20px; border-radius: 5px; margin: 20px 0; }}
        .metric {{ display: inline-block; margin: 10px 20px; }}
        .metric-label {{ font-weight: bold; color: #666; }}
        .metric-value {{ font-size: 24px; color: #2196F3; }}
        table {{ width: 100%; border-collapse: collapse; margin: 20px 0; }}
        th {{ background: #4CAF50; color: white; padding: 12px; text-align: left; }}
        td {{ padding: 10px; border-bottom: 1px solid #ddd; }}
        tr:hover {{ background: #f5f5f5; }}
        .good {{ color: #4CAF50; font-weight: bold; }}
        .warning {{ color: #FF9800; font-weight: bold; }}
        .bad {{ color: #F44336; font-weight: bold; }}
        .diagnostics {{ background: #fff3e0; padding: 15px; border-left: 4px solid #FF9800; margin: 20px 0; }}
    </style>
</head>
<body>
    <div class="container">
        <h1>ISMRM 2015 Tractography Validation Report</h1>
        <p><strong>Generated:</strong> {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</p>
        <p><strong>NIM File:</strong> {nim_file}</p>
        <p><strong>Tracks File:</strong> {tracks_file}</p>

        <div class="summary">
            <h2>Summary Statistics</h2>
            <div class="metric">
                <div class="metric-label">Valid Bundles</div>
                <div class="metric-value">{valid_bundles} / {total_bundles}</div>
            </div>
            <div class="metric">
                <div class="metric-label">Total Streamlines</div>
                <div class="metric-value">{total_user_tracks:,}</div>
            </div>
            <div class="metric">
                <div class="metric-label">Average Coverage</div>
                <div class="metric-value">{avg_coverage:.1f}%</div>
            </div>
            <div class="metric">
                <div class="metric-label">Average Overreach</div>
                <div class="metric-value">{avg_overreach:.1f}%</div>
            </div>
        </div>

        <h2>Per-Bundle Results</h2>
        <table>
            <tr>
                <th>Bundle</th>
                <th>User Tracks</th>
                <th>GT Tracks</th>
                <th>Coverage (%)</th>
                <th>Overreach (%)</th>
                <th>Avg Length (mm)</th>
                <th>Status</th>
            </tr>
"""

        for bundle_name in sorted(self.results['bundles'].keys()):
            metrics = self.results['bundles'][bundle_name]

            status_class = 'good' if metrics['valid'] and metrics['coverage'] > 50 else 'warning' if metrics['valid'] else 'bad'
            status_text = '✓ Valid' if metrics['valid'] else '✗ Empty'

            html += f"""
            <tr>
                <td><strong>{bundle_name}</strong></td>
                <td>{metrics['user_count']}</td>
                <td>{metrics['gt_count']}</td>
                <td>{metrics['coverage']:.1f}%</td>
                <td>{metrics['overreach']:.1f}%</td>
                <td>{metrics['avg_length_user']:.1f}</td>
                <td class="{status_class}">{status_text}</td>
            </tr>
"""

        html += """
        </table>

        <div class="diagnostics">
            <h2>Coordinate System Diagnostics</h2>
"""

        if 'input_coords' in self.results['diagnostics']:
            ic = self.results['diagnostics']['input_coords']
            html += f"""
            <p><strong>Input Voxel Coordinates:</strong></p>
            <ul>
                <li>X: [{ic['x_range'][0]:.1f}, {ic['x_range'][1]:.1f}]</li>
                <li>Y: [{ic['y_range'][0]:.1f}, {ic['y_range'][1]:.1f}]</li>
                <li>Z: [{ic['z_range'][0]:.1f}, {ic['z_range'][1]:.1f}]</li>
            </ul>
"""

        if 'world_coords' in self.results['diagnostics']:
            wc = self.results['diagnostics']['world_coords']
            html += f"""
            <p><strong>World Coordinates (after transform):</strong></p>
            <ul>
                <li>X: [{wc['x_range'][0]:.1f}, {wc['x_range'][1]:.1f}] mm</li>
                <li>Y: [{wc['y_range'][0]:.1f}, {wc['y_range'][1]:.1f}] mm</li>
                <li>Z: [{wc['z_range'][0]:.1f}, {wc['z_range'][1]:.1f}] mm</li>
            </ul>
"""

        html += """
        </div>

        <h2>Interpretation</h2>
        <ul>
            <li><strong>Coverage:</strong> Percentage of ground truth streamlines covered by user tractography (higher is better)</li>
            <li><strong>Overreach:</strong> Percentage of user streamlines not matching ground truth (lower is better)</li>
            <li><strong>Valid Bundles:</strong> Bundles with at least one streamline passing ROI constraints</li>
            <li><strong>Good Performance:</strong> Coverage > 50%, Overreach < 30%</li>
        </ul>

    </div>
</body>
</html>
"""

        with open(report_path, 'w') as f:
            f.write(html)

        print(f"✓ Report saved: {report_path}")

        # Save JSON results
        json_path = self.output_dir / 'validation_results.json'
        with open(json_path, 'w') as f:
            json.dump(self.results, f, indent=2)
        print(f"✓ Results saved: {json_path}")

        return report_path

    def validate(self, nim_file, tracks_file):
        """
        Main validation workflow

        Steps:
        1. Load HINEC tracks (voxel space)
        2. Segment bundles using ROI constraints (voxel space)
        3. Convert to world space for distance metrics
        4. Compare with ground truth
        5. Compute metrics
        6. Generate report
        """
        print("=" * 60)
        print("ISMRM TRACTOGRAPHY VALIDATION")
        print("=" * 60)

        # Load tracks (in voxel space)
        voxel_tracks = self.load_hinec_tracks(tracks_file)

        # Get reference for later world space conversion
        reference_nii = self.scoring_dir / 't1.nii.gz'
        if not reference_nii.exists():
            raise FileNotFoundError(f"Reference T1 not found: {reference_nii}")

        # Segment and evaluate each bundle
        print(f"\n=== Bundle Segmentation & Evaluation ===")
        for bundle_name, bundle_config in self.seg_config.items():
            # Segment user tractography (ROI filtering in voxel space)
            user_bundle_voxel = self.segment_bundle_roi(voxel_tracks, bundle_name, bundle_config)

            # Convert segmented bundle to world space for distance metrics
            if len(user_bundle_voxel) > 0:
                user_bundle_world = self.convert_to_world_space(user_bundle_voxel, reference_nii)
            else:
                user_bundle_world = []

            # Load ground truth
            gt_bundle = self.load_ground_truth(bundle_name)

            # Compute metrics (uses world space for distance calculations)
            metrics = self.compute_bundle_metrics(user_bundle_world, gt_bundle, bundle_name)
            self.results['bundles'][bundle_name] = metrics

            # Save segmented bundle
            if len(user_bundle_world) > 0:
                bundle_path = self.output_dir / f'{bundle_name}_segmented.trk'
                tractogram = Tractogram(streamlines=user_bundle_world)
                try:
                    save_trk(tractogram, str(bundle_path))
                except Exception as e:
                    print(f"      Warning: Could not save bundle TRK: {e}")

        # Generate report
        report_path = self.generate_report(nim_file, tracks_file)

        print(f"\n{'=' * 60}")
        print("VALIDATION COMPLETE")
        print(f"{'=' * 60}")
        print(f"Report: {report_path}")
        print(f"\nSummary:")
        print(f"  Valid bundles: {self.results['summary']['valid_bundles']} / {self.results['summary']['total_bundles']}")
        print(f"  Average coverage: {self.results['summary']['avg_coverage']:.1f}%")
        print(f"  Average overreach: {self.results['summary']['avg_overreach']:.1f}%")

        return self.results


def main():
    parser = argparse.ArgumentParser(
        description="Validate HINEC tractography against ISMRM ground truth",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.add_argument('--nim-file', required=True,
                       help='HINEC nim structure (.mat file)')
    parser.add_argument('--tracks-file', required=True,
                       help='HINEC tractography results (.mat file)')
    parser.add_argument('--scoring-dir', required=True,
                       help='ISMRM scoring data directory')
    parser.add_argument('--output-dir', default='validation_results',
                       help='Output directory for results (default: validation_results)')

    args = parser.parse_args()

    try:
        validator = ISMRMValidator(args.scoring_dir, args.output_dir)
        results = validator.validate(args.nim_file, args.tracks_file)

        sys.exit(0)

    except Exception as e:
        print(f"\nError: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == '__main__':
    main()
