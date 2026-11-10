"""Minimal YAML reader mirroring src/nim_utils/load_config_yaml.m.

The MATLAB parser only understands two-level `section: / key: value` files with
`#` comments, so we mirror exactly that (no PyYAML dependency) and apply the same
defaults + validation as validate_config()/validate_tractography_params().
"""
from __future__ import annotations

PREPROCESSING_DEFAULTS = dict(
    run_denoising=True, denoise_method="dwidenoise", run_motion_correction=True,
    run_eddy=True, improve_mask=True, atlas_type="jhu", t1_available=False,
    use_t1_registration=False, register_to_mni=False,
)

TRACTOGRAPHY_DEFAULTS = dict(
    algorithm="hinec", integration_order=4, interp_method="trilinear", step_size=0.2,
    adaptive_step=False, rkf_tolerance=0.01, rkf_safety=0.9, step_min=0.01, step_max=1.0,
    seed_density=4, seed_strategy="uniform", termination_fa=0.15, angle_thresh=35,
    max_steps=1000, min_length=35, act_enabled=True, enable_diagnostics=True,
)


def _parse_value(s: str):
    s = s.strip()
    if (s.startswith("'") and s.endswith("'")) or (s.startswith('"') and s.endswith('"')):
        return s[1:-1]
    low = s.lower()
    if low == "true":
        return True
    if low == "false":
        return False
    if low in ("null", "~"):
        return None
    try:
        v = float(s)
        return int(v) if v.is_integer() and "." not in s and "e" not in low else v
    except ValueError:
        return s


def load_config_yaml(path: str) -> dict:
    cfg: dict = {}
    section = None
    with open(path) as fh:
        for raw in fh:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            if line.endswith(":") and " " not in line:
                section = line[:-1]
                cfg[section] = {}
                continue
            if ":" in line:
                key, _, val = line.partition(":")
                val = val.split("#", 1)[0].strip()
                if not val:
                    continue
                target = cfg if section is None else cfg[section]
                target[key.strip()] = _parse_value(val)

    cfg.setdefault("preprocessing", {})
    cfg.setdefault("tractography", {})
    for k, v in PREPROCESSING_DEFAULTS.items():
        cfg["preprocessing"].setdefault(k, v)
    for k, v in TRACTOGRAPHY_DEFAULTS.items():
        cfg["tractography"].setdefault(k, v)
    _validate(cfg["tractography"])
    return cfg


def _validate(p: dict) -> None:
    if p["step_size"] <= 0:
        raise ValueError(f"step_size must be positive, got {p['step_size']}")
    if not (0 < p["angle_thresh"] <= 180):
        raise ValueError(f"angle_thresh must be in (0,180], got {p['angle_thresh']}")
    if p["integration_order"] not in (1, 2, 4, 5):
        raise ValueError(f"integration_order must be 1,2,4,5, got {p['integration_order']}")
    if p["integration_order"] == 5:
        if p["rkf_tolerance"] <= 0:
            raise ValueError("rkf_tolerance must be positive")
        if p["step_min"] <= 0 or p["step_min"] >= p["step_max"]:
            raise ValueError("invalid step bounds")
        if not (0 < p["rkf_safety"] <= 1):
            raise ValueError("rkf_safety must be in (0,1]")
    if not (0 <= p["termination_fa"] <= 1):
        raise ValueError("termination_fa must be in [0,1]")
