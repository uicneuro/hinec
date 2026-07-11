#!/usr/bin/env python
"""run_pft_dipy.py — Particle Filtering Tractography (Girard 2014) via DIPY.

The winning ISMRM-2015 approach per Renauld 2023: CSD FOD + PROBABILISTIC tracking
+ WM/GM/CSF anatomical priors (CMC) + particle-filter BACKTRACKING (when a
streamline stops somewhere invalid, back up and find an alternative until it
terminates in gray matter). Targets the recall/overlap deficit of deterministic
tracking. Outputs a .trk (RAS world) for the same scilpy Renauld-2023 scorer.

CSD is computed in DIPY's native SH basis so the direction getter is correct
(avoids the gradient/basis pitfalls that broke the MRtrix ACT runs). Tissue masks
are DWI-space (from HINEC preprocessing), so no T1->DWI registration to misalign.

Usage:
  run_pft_dipy.py --out tracks.trk [--density 1] [--seed wm|interface]
                  [--rng K]  (K fixes numpy RNG for reproducibility; omit = stochastic)
"""
import argparse, sys, time
import numpy as np
import nibabel as nib
from dipy.io.gradients import read_bvals_bvecs
from dipy.core.gradients import gradient_table
from dipy.reconst.csdeconv import ConstrainedSphericalDeconvModel, auto_response_ssst
from dipy.direction import ProbabilisticDirectionGetter
from dipy.data import default_sphere
from dipy.tracking.stopping_criterion import CmcStoppingCriterion
from dipy.tracking.local_tracking import ParticleFilteringTracking
from dipy.tracking.streamline import Streamlines
from dipy.tracking import utils
from dipy.io.stateful_tractogram import StatefulTractogram, Space
from dipy.io.streamline import save_trk

DATA = "data/ismrm2015"

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--density", type=int, default=1, help="seeds per voxel per axis")
    ap.add_argument("--seed", choices=["wm", "interface"], default="interface")
    ap.add_argument("--sh_order", type=int, default=6)
    ap.add_argument("--max_angle", type=float, default=20.0)
    ap.add_argument("--step", type=float, default=0.5, help="mm")
    ap.add_argument("--minlen", type=float, default=20.0)
    ap.add_argument("--maxlen", type=float, default=250.0)
    ap.add_argument("--particles", type=int, default=15)
    ap.add_argument("--rng", type=int, default=-1, help=">=0 fixes numpy seed (reproducibility)")
    args = ap.parse_args()
    if args.rng >= 0:
        np.random.seed(args.rng)

    t0 = time.time()
    dwi = nib.load(f"{DATA}/ismrm2015.nii.gz")
    data = dwi.get_fdata(); affine = dwi.affine; vox = float(dwi.header.get_zooms()[0])
    bvals, bvecs = read_bvals_bvecs(f"{DATA}/ismrm2015.bval", f"{DATA}/ismrm2015.bvec")
    gtab = gradient_table(bvals, bvecs)
    # Proper T1-derived tissue PVE (FSL FAST on the anatomical T1, DWI space).
    # Falls back to the old FA-tertile masks only if the PVE maps are absent.
    import os
    if os.path.exists(f"{DATA}/ismrm2015_WM_pve.nii.gz"):
        wm = nib.load(f"{DATA}/ismrm2015_WM_pve.nii.gz").get_fdata()
        gm = nib.load(f"{DATA}/ismrm2015_GM_pve.nii.gz").get_fdata()
        csf = nib.load(f"{DATA}/ismrm2015_CSF_pve.nii.gz").get_fdata()
        tot = wm + gm + csf; m = tot > 1e-3
        wm[m] /= tot[m]; gm[m] /= tot[m]; csf[m] /= tot[m]
    else:
        wm = (nib.load(f"{DATA}/ismrm2015_WM_mask.nii.gz").get_fdata() > 0.5).astype(float)
        gm = (nib.load(f"{DATA}/ismrm2015_GM_mask.nii.gz").get_fdata() > 0.5).astype(float)
        csf = (nib.load(f"{DATA}/ismrm2015_CSF_mask.nii.gz").get_fdata() > 0.5).astype(float)
    brain = (wm + gm + csf) > 0.1
    print(f"[{time.time()-t0:.0f}s] loaded DWI {data.shape}, tissue WM={int((wm>0.5).sum())} GM={int((gm>0.5).sum())} CSF={int((csf>0.5).sum())}", flush=True)

    # --- CSD in DIPY's native basis (correct direction getter) ---
    resp, ratio = auto_response_ssst(gtab, data, roi_radii=10, fa_thr=0.7)
    print(f"[{time.time()-t0:.0f}s] CSD response {resp[0]} ratio {ratio:.3f}; fitting FOD...", flush=True)
    try:
        csd = ConstrainedSphericalDeconvModel(gtab, resp, sh_order_max=args.sh_order)
    except TypeError:
        csd = ConstrainedSphericalDeconvModel(gtab, resp, sh_order=args.sh_order)
    fit = csd.fit(data, mask=brain)
    shm = fit.shm_coeff
    print(f"[{time.time()-t0:.0f}s] FOD fit done; SH shape {shm.shape}", flush=True)

    dg = ProbabilisticDirectionGetter.from_shcoeff(shm, max_angle=args.max_angle, sphere=default_sphere)

    # --- CMC anatomical stopping criterion from tissue PVE (binary here) ---
    sc = CmcStoppingCriterion.from_pve(wm, gm, csf, step_size=args.step, average_voxel_size=vox)

    # --- seeds: GM-WM interface (PFT standard) or WM ---
    if args.seed == "interface":
        from scipy.ndimage import binary_dilation
        seed_mask = (binary_dilation(gm > 0.5) & (wm > 0.5))
    else:
        seed_mask = wm > 0.5
    seeds = utils.seeds_from_mask(seed_mask, affine, density=args.density)
    print(f"[{time.time()-t0:.0f}s] {seeds.shape[0]} seeds ({args.seed}); tracking PFT...", flush=True)

    pft = ParticleFilteringTracking(
        dg, sc, seeds, affine, step_size=args.step,
        max_cross=1, maxlen=int(args.maxlen/args.step),
        pft_back_tracking_dist=2, pft_front_tracking_dist=1,
        particle_count=args.particles, return_all=False)
    sl = Streamlines(pft)
    # length filter
    from dipy.tracking.utils import length as sl_length
    lens = np.array(list(sl_length(sl)))
    keep = (lens >= args.minlen) & (lens <= args.maxlen)
    sl = Streamlines([s for s, k in zip(sl, keep) if k])
    print(f"[{time.time()-t0:.0f}s] {len(sl)} streamlines (after {args.minlen}-{args.maxlen}mm filter)", flush=True)

    # save with the SCORING T1 as reference (streamlines are already RAS world) so
    # the scorer's ROIs (T1 space) are space-compatible.
    t1 = nib.load(f"{DATA}/scoring_data_Renauld2023/t1.nii.gz")
    sft = StatefulTractogram(sl, t1, Space.RASMM)
    save_trk(sft, args.out, bbox_valid_check=False)
    print(f"[{time.time()-t0:.0f}s] saved {args.out} ({len(sl)} streamlines)", flush=True)

if __name__ == "__main__":
    main()
