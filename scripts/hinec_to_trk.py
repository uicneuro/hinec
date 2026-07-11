#!/usr/bin/env python3
"""
hinec_to_trk.py - Convert HINEC tractography to TrackVis TRK format

Converts MATLAB .mat files with HINEC tracks to TRK format compatible with
ISMRM 2015 scoring scripts and Dipy StatefulTractogram.

CRITICAL: Handles coordinate system conversion and proper header setup to avoid
the 0.5 voxel shift issue mentioned in ISMRM scoring documentation.

Usage:
    python hinec_to_trk.py tracks.mat reference.nii.gz output.trk [--space voxmm]

Arguments:
    tracks_mat      - HINEC tractography results (.mat file)
    reference_nii   - Reference NIfTI image (for header/affine)
    output_trk      - Output TrackVis file
    --space         - Coordinate space: 'vox', 'voxmm', or 'rasmm' (default: voxmm)
"""

import argparse
import sys
import numpy as np

try:
    import nibabel as nib
    from nibabel.streamlines import Tractogram, save
except ImportError:
    print("Error: nibabel not installed")
    print("Install with: pip install nibabel")
    sys.exit(1)

try:
    import scipy.io
except ImportError:
    print("Error: scipy not installed")
    print("Install with: pip install scipy")
    sys.exit(1)

try:
    import h5py
    H5PY_AVAILABLE = True
except ImportError:
    H5PY_AVAILABLE = False


def _load_hdf5_tracks(mat_file):
    """Load tracks from MATLAB v7.3 (HDF5) file via h5py."""
    with h5py.File(mat_file, 'r') as f:
        if 'tracks' in f:
            tracks_ref = f['tracks']
        elif 'tracts' in f:
            tracks_ref = f['tracts']
        else:
            raise ValueError("No 'tracks' or 'tracts' field in HDF5 file")

        tracks = []
        if tracks_ref.dtype == np.dtype('O'):
            num_tracks = tracks_ref.shape[0] * tracks_ref.shape[1]
            tracks_flat = tracks_ref[...].flatten()
            for i in range(num_tracks):
                track = f[tracks_flat[i]][()]
                if track.size == 0:
                    continue
                if track.shape[0] == 3 and track.shape[1] != 3:
                    track = track.T  # 3xN -> Nx3
                elif track.shape[1] != 3:
                    continue
                tracks.append(track.astype(np.float32) - 1.0)
        else:
            data = tracks_ref[()]
            if data.ndim == 2 and data.shape[1] == 3:
                tracks.append(data.astype(np.float32) - 1.0)
            else:
                raise ValueError(f"Unexpected track data shape: {data.shape}")
    return tracks


def load_hinec_tracks(mat_file):
    """
    Load HINEC tracks from MATLAB .mat file (handles both standard and v7.3/HDF5).

    Expected structure:
        tracks: cell array where each cell is an Nx3 matrix
        Each track is stored as [X, Y, Z] coordinates in voxel space (1-based MATLAB)

    Returns:
        list of numpy arrays (Nx3 each)
    """
    print(f"Loading HINEC tracks from {mat_file}...")

    # MATLAB v7.3 is HDF5, which scipy.io cannot read. Depending on the scipy
    # version this surfaces as NotImplementedError OR MatReadError ("first 20
    # bytes == 0") -- so fall back to the h5py loader on ANY scipy read failure
    # if the file is in fact a valid HDF5/v7.3 file.
    try:
        mat_data = scipy.io.loadmat(mat_file)
    except Exception as e:  # noqa: BLE001 -- broad on purpose (scipy raises several types)
        is_hdf5 = H5PY_AVAILABLE and h5py.is_hdf5(mat_file)
        if is_hdf5 or 'v7.3' in str(e) or 'HDF' in str(e):
            if not H5PY_AVAILABLE:
                raise ImportError(
                    "This is a MATLAB v7.3 file (HDF5). Install h5py: pip install h5py"
                )
            print("  Format: MATLAB v7.3 (HDF5)")
            tracks = _load_hdf5_tracks(mat_file)
            print(f"  Loaded {len(tracks)} tracks")
            if len(tracks) > 0:
                all_coords = np.vstack(tracks)
                print(f"  Coordinate range (0-based voxels):")
                print(f"    X: [{all_coords[:,0].min():.1f}, {all_coords[:,0].max():.1f}]")
                print(f"    Y: [{all_coords[:,1].min():.1f}, {all_coords[:,1].max():.1f}]")
                print(f"    Z: [{all_coords[:,2].min():.1f}, {all_coords[:,2].max():.1f}]")
            return tracks
        raise

    # Standard MATLAB format
    if 'tracks' in mat_data:
        tracks_cell = mat_data['tracks']
    elif 'tracts' in mat_data:
        tracks_cell = mat_data['tracts']
    else:
        raise ValueError("No 'tracks' or 'tracts' field found in .mat file")

    tracks = []
    for i in range(tracks_cell.shape[0]):
        track = tracks_cell[i, 0]
        if track.size > 0:
            track_0based = track.astype(np.float32) - 1.0
            tracks.append(track_0based)

    print(f"  Loaded {len(tracks)} tracks")

    if len(tracks) > 0:
        # Show coordinate statistics
        all_coords = np.vstack(tracks)
        print(f"  Coordinate range (0-based voxels):")
        print(f"    X: [{all_coords[:,0].min():.1f}, {all_coords[:,0].max():.1f}]")
        print(f"    Y: [{all_coords[:,1].min():.1f}, {all_coords[:,1].max():.1f}]")
        print(f"    Z: [{all_coords[:,2].min():.1f}, {all_coords[:,2].max():.1f}]")

    return tracks


def convert_to_trk(tracks, reference_nii, output_trk, space='voxmm', trk_ref=None):
    """
    Convert HINEC tracks to TRK format with proper header.

    Args:
        tracks: List of Nx3 numpy arrays (0-based voxel coordinates)
        reference_nii: Path to NIfTI whose affine defines the voxel grid the
            tracks came from (used for voxel -> world conversion).
        output_trk: Output TRK file path
        space: Coordinate space ('vox', 'voxmm', or 'rasmm')
        trk_ref: Optional separate NIfTI to attach as the saved TRK's reference.
            Use this when the ROI/atlas masks the scorer compares against live
            in a different reference space (e.g. scoring T1 at 1 mm) than the
            DWI the tracks came from (e.g. 2 mm). World mm coordinates remain
            physically correct; only the saved reference metadata changes.

    Coordinate spaces:
        'vox'    - Voxel indices (0-based)
        'voxmm'  - Voxel corner coordinates in mm
        'rasmm'  - RAS+ world coordinates in mm
    """
    print(f"\nLoading reference NIfTI: {reference_nii}")
    ref_img = nib.load(reference_nii)
    affine = ref_img.affine

    print(f"  Dimensions: {ref_img.shape}")
    print(f"  Voxel size: {ref_img.header.get_zooms()[:3]} mm")
    print(f"  Affine transform:")
    print(f"    {affine}")

    # Convert tracks to requested space
    print(f"\nConverting tracks to '{space}' space...")

    if space == 'vox':
        # Keep as 0-based voxel indices
        converted_tracks = tracks
        print("  Using voxel indices (0-based)")

    elif space == 'voxmm':
        # Convert to voxel corner coordinates (voxel indices * voxel_size)
        # This matches TrackVis convention
        voxel_sizes = np.array(ref_img.header.get_zooms()[:3])
        converted_tracks = [track * voxel_sizes for track in tracks]
        print(f"  Scaled by voxel sizes: {voxel_sizes}")

    elif space == 'rasmm':
        # Convert to RAS+ world coordinates using affine
        converted_tracks = []
        for track in tracks:
            # Add homogeneous coordinate
            track_homog = np.column_stack([track, np.ones(len(track))])
            # Apply affine transform
            track_world = track_homog @ affine.T
            # Remove homogeneous coordinate
            converted_tracks.append(track_world[:, :3].astype(np.float32))
        print("  Transformed to RAS+ world coordinates")

    else:
        raise ValueError(f"Unknown space: {space}. Use 'vox', 'voxmm', or 'rasmm'")

    # Show converted coordinate statistics
    if len(converted_tracks) > 0:
        all_coords = np.vstack(converted_tracks)
        print(f"  Converted coordinate range:")
        print(f"    X: [{all_coords[:,0].min():.1f}, {all_coords[:,0].max():.1f}]")
        print(f"    Y: [{all_coords[:,1].min():.1f}, {all_coords[:,1].max():.1f}]")
        print(f"    Z: [{all_coords[:,2].min():.1f}, {all_coords[:,2].max():.1f}]")

    # Save TRK. If a separate trk_ref is provided, use dipy's StatefulTractogram
    # so the saved TRK header carries the trk_ref's affine/dims/voxel sizes —
    # what scilpy compares against the ROI masks. The streamlines are already
    # in physical mm world coords (computed using reference_nii's affine), so
    # they remain at the correct anatomical location regardless of which
    # reference is attached for metadata purposes.
    if space == 'rasmm' and trk_ref is not None:
        try:
            from dipy.io.stateful_tractogram import StatefulTractogram, Space
            from dipy.io.streamline import save_tractogram
        except ImportError:
            print("Error: dipy is required for --trk-ref. Install with: pip install dipy")
            sys.exit(1)
        print(f"\nAttaching saved TRK reference: {trk_ref}")
        trk_ref_img = nib.load(trk_ref)
        print(f"  TRK header dims: {trk_ref_img.shape[:3]}, voxel size: {trk_ref_img.header.get_zooms()[:3]} mm")
        sft = StatefulTractogram(converted_tracks, trk_ref, Space.RASMM)
        print(f"Saving to {output_trk}...")
        save_tractogram(sft, output_trk, bbox_valid_check=False)
    else:
        print(f"\nCreating Tractogram...")
        tractogram = Tractogram(
            streamlines=converted_tracks,
            affine_to_rasmm=affine if space != 'rasmm' else np.eye(4)
        )
        print(f"Saving to {output_trk}...")
        save(tractogram, output_trk)

    print("✓ Conversion complete!")
    print(f"\nTo verify alignment:")
    print(f"  1. Load in TrackVis or MI-Brain")
    print(f"  2. Overlay with {reference_nii}")
    print(f"  3. Check if tracks follow anatomical structures")


def main():
    parser = argparse.ArgumentParser(
        description="Convert HINEC tractography to TrackVis TRK format",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Convert to voxel-mm space (recommended for ISMRM scoring)
  python hinec_to_trk.py tracks.mat t1.nii.gz output.trk --space voxmm

  # Convert to RAS+ world space
  python hinec_to_trk.py tracks.mat t1.nii.gz output.trk --space rasmm

  # Keep as voxel indices
  python hinec_to_trk.py tracks.mat t1.nii.gz output.trk --space vox

Coordinate spaces:
  vox    - Voxel indices (0-based, matches HINEC internal format)
  voxmm  - Voxel corner coordinates in mm (TrackVis convention)
  rasmm  - RAS+ world coordinates in mm (NIFTI convention)

Note: ISMRM scoring expects proper coordinate space attributes in TRK header.
      The default 'voxmm' should work with Dipy's StatefulTractogram.
        """
    )

    parser.add_argument('tracks_mat', help='HINEC tractography .mat file')
    parser.add_argument('reference_nii', help='Reference NIfTI image')
    parser.add_argument('output_trk', help='Output TrackVis .trk file')
    parser.add_argument('--space', default='rasmm',
                       choices=['vox', 'voxmm', 'rasmm'],
                       help='Coordinate space (default: rasmm — what the ISMRM 2023 scilpy scorer expects)')
    parser.add_argument('--trk-ref', default=None,
                       help='Optional separate NIfTI to attach as the saved TRK reference. Use when '
                            'the scorer compares against masks in a different reference (e.g. scoring '
                            'T1 at 1 mm) than the DWI where tracks live (e.g. 2 mm). Only used with '
                            '--space rasmm.')

    args = parser.parse_args()

    try:
        # Load HINEC tracks
        tracks = load_hinec_tracks(args.tracks_mat)

        if len(tracks) == 0:
            print("Warning: No tracks found in .mat file!")
            sys.exit(1)

        # Convert to TRK
        convert_to_trk(tracks, args.reference_nii, args.output_trk, args.space,
                       trk_ref=args.trk_ref)

    except Exception as e:
        print(f"\nError: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == '__main__':
    main()
