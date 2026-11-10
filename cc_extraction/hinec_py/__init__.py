"""hinec_py — a faithful NumPy/SciPy port of the HINEC MATLAB pipeline.

Modules
-------
dti       : nim_read → brain mask → nim_dt_spd → nim_eig → nim_fa → mask improvement → tissue masks
tracking  : nim_tractography_hinec (lattice seeding, Euler/RK2/RK4/RKF45, trilinear/cubic interpolation, ACT)
cc        : corpus-callosum ROI construction and streamline selection
io        : NIfTI / TRK / MATLAB (.mat) writers compatible with the HINEC viewers
config    : minimal reader for HINEC's YAML presets (same subset the MATLAB parser understands)
"""
__version__ = "0.1.0"
