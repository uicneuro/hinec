# RKF45 Adaptive Stepping — Usage

`integrator.method: rkf45` runs the HINEC tracker with a Dormand–Prince embedded
Runge–Kutta pair and, by default, adaptive step-size control: the step shrinks where the
direction field turns sharply and grows where it is straight, holding the estimated local
error at a tolerance instead of at a step size you guessed. The scheme itself is derived in
[RKF45](RKF.md).

The same machinery is available in the MMF tracker (`integrator.method: rkf45` with
`algorithm: mmf`), applied to the coupled frame system rather than to a position alone.

---

## Enabling it

```yaml
tractography:
  integrator:
    method: rkf45
    step: 0.2        # INITIAL step, in voxels
    tolerance: 0.02  # accepted local error, in voxels
    step_min: 0.02
    step_max: 0.5
```

`config/hinec_dti.yml` and `config/hinec_csd.yml` ship with exactly these integrator
settings; `config/hinec_dti_cubic.yml` selects `rkf45` and leaves the rest at the schema
defaults. To sweep one knob without writing a config:

```bash
./bin/run_tractography.sh hinec_dti --score --set integrator.tolerance=0.005
```

---

## Parameters

All live under `tractography.integrator`; defaults are from
`src/nim_utils/nim_config_schema.m`.

| Key | Default | Applies to | Meaning |
|---|---|---|---|
| `method` | `rk4` | all | `euler` \| `rk2` \| `rk4` \| `rkf45` |
| `step` | `0.2` | all | fixed step, or the **initial** step for `rkf45`, in voxels |
| `adaptive` | `true` | `rkf45` only | step-size control; `false` gives fixed-step RKF45 |
| `tolerance` | `0.01` | `rkf45` only | accepted local error per step, in voxels |
| `step_min` | `0.01` | `rkf45` only | lower clamp on the adapted step |
| `step_max` | `1.0` | `rkf45` only | upper clamp on the adapted step |
| `safety` | `0.9` | `rkf45` only | safety factor in the step-size update, \(< 1\) |

Settings marked `rkf45` only are ignored by the fixed-step methods, by design — leaving them
in a config that switches to `rk4` is harmless.

### How the knobs interact

- **`tolerance`** is a *local* error bound per step, in voxels, not a bound on the final
  track. Tightening it reduces the step and raises the number of steps and interpolations;
  loosening it does the reverse. It is the main accuracy/cost dial.
- **`step`** only sets where adaptation starts. Too large an initial step spends the first
  few steps being rejected; too small merely costs a few extra accepted steps.
- **`step_min` / `step_max`** clamp the adapted step. A `step_min` that binds turns off error
  control in exactly the regions that needed it, so if runs are slow, prefer loosening
  `tolerance` over raising `step_min`.
- Step growth is additionally capped at 2× per step, so recovery from a rejected step is
  gradual by construction.
- **`adaptive: false`** keeps the seven-stage pair but fixes the step. This is a real third
  mode, kept so that the tableau can be compared independently of the step control.

---

## Cost

Each attempted step costs seven interpolations of the direction field against four for RK4,
and rejected steps are attempted work that is thrown away. Against that, adaptation removes
steps in straight regions. The net cost per track is therefore data-dependent; measure it on
your own volume rather than assuming a fixed ratio.

Accuracy, on the other hand, is not decided by the tableau alone. The order any integrator
can realise here is capped by the smoothness of the interpolated field, so pairing a
high-order scheme with `interpolation.method: trilinear` (\(C^0\)) buys much less than the
stage count suggests. See [Solution Verification](CONVERGENCE.md).

---

## Diagnostics

With `tractography.diagnostics: true` (the default) a run reports, alongside the usual
timing breakdown:

- the integration method and whether adaptive stepping is enabled, with the tolerance;
- total integration steps, average steps per track, and steps per second;
- step **rejection** and **retry** counts — but only when they are being accumulated: the
  parallel tracking path does not aggregate them across workers, so the line is absent from
  a normal `parfor` run;
- the termination analysis — how many track ends fell to the FA threshold or the angle
  threshold, how many had no initial direction, and, with ACT active, the GM / CSF / outside
  counts.

The termination analysis is the quickest way to tell whether the integrator is what limits a
run. If most ends fall to `FA threshold` or `Angle threshold`, tightening the tolerance will
not change the result.

---

## Troubleshooting

| Symptom | Likely cause | First thing to try |
|---|---|---|
| Many rejected steps | tolerance tight relative to the initial step | raise `integrator.tolerance`, or lower `integrator.step` |
| Very slow runs | tolerance tighter than the data supports | raise `integrator.tolerance`; check whether `step_min` is binding |
| Tracks end too early | a termination criterion, not the integrator | check `termination.max_arc` (arc length in voxels — `max_steps` is derived from it), `termination.fa_min`, and `termination.angle_max` |
| Tracks hairpin near their ends | turn budget inert at this step | `termination.angle_max × integrator.step` must stay below 90°, or the criterion cannot fire |

---

## Comparing against fixed-step RK4

Run the same config twice, changing only the integrator, and compare in a scorable run dir
rather than by eye:

```bash
./bin/run_tractography.sh hinec_dti --score --set integrator.method=rk4
./bin/run_tractography.sh hinec_dti --score --set integrator.method=rkf45
```

Each gets its own tagged run directory; compare `scoring/renauld2023/results.json`. For a
controlled convergence study rather than a scoring comparison, follow
[Solution Verification](CONVERGENCE.md), which pins the seeds and disables the angle
criterion so that only the integrator varies.

---

## References

1. **Dormand, J. R., & Prince, P. J.** (1980). A family of embedded Runge-Kutta formulae.
   *Journal of Computational and Applied Mathematics*, 6(1), 19–26.
2. **Hairer, E., Nørsett, S. P., & Wanner, G.** (1993). *Solving Ordinary Differential
   Equations I: Nonstiff Problems* (2nd ed.). Springer.
3. **Basser, P. J., Pajevic, S., Pierpaoli, C., Duda, J., & Aldroubi, A.** (2000). In vivo
   fiber tractography using DT-MRI data. *Magnetic Resonance in Medicine*, 44(4), 625–632.
