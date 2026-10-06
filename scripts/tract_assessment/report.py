"""Write an assessment's outputs: summary.json, pairs.csv, figures and a short
HTML report that leads with the reading, not the numbers."""
import base64
import csv
import html
import json
import os
from pathlib import Path

import numpy as np

# Reading bands. These are CONVENTIONS for turning numbers into words, not
# validated pass/fail thresholds (the tract-geometry study found no universal
# cutoff). They are printed in every report so a reader can disagree with them.
BANDS = {
    'overlap':   (0.8, 0.5),     # Jaccard / Dice: >=0.8 high, >=0.5 moderate, else low
    'rank':      (0.8, 0.5),     # Spearman
    'shift':     (0.02, 0.10),   # relative median change: <=2% high, <=10% moderate
}


def band(kind, value):
    if value is None or not np.isfinite(value):
        return 'n/a'
    hi, mid = BANDS[kind]
    if kind == 'shift':
        return 'high' if value <= hi else 'moderate' if value <= mid else 'low'
    return 'high' if value >= hi else 'moderate' if value >= mid else 'low'


def reading(h):
    """Plain-language reading of one headline dict, tier by tier."""
    spatial = band('overlap', h['occupied_jaccard'])
    endpoints = band('overlap', h['endpoint_jaccard'])
    geom_shift = band('shift', h['kappa_w_rel_median_change'])
    geom_rank = band('rank', h['kappa_w_paired_spearman'])
    lines = [
        f"Spatial agreement is <b>{spatial}</b> (occupied-voxel Jaccard {_f(h['occupied_jaccard'])}); "
        f"endpoint agreement is <b>{endpoints}</b> ({_f(h['endpoint_jaccard'])}).",
        f"Matched streamlines ({h['pairs']} pairs, matched by {h['correspondence']}) are separated by a median "
        f"{_f(h['path_separation_median_mm'])} mm; the worst 5% by {_f(h['path_separation_tail_mm'])} mm or more.",
        f"MMF curvature κ<sub>w</sub>: the bundle-level distribution agreement is <b>{geom_shift}</b> "
        f"(median shift {_pct(h['kappa_w_rel_median_change'])}); per-streamline agreement is <b>{geom_rank}</b> "
        f"(paired Spearman {_f(h['kappa_w_paired_spearman'])}).",
    ]
    if geom_shift == 'high' and spatial != 'high':
        lines.append('<b>Warning:</b> geometry agrees but location does not. Similar κ<sub>w</sub> statistics do not '
                     'show that the tractograms agree; report the spatial result alongside them.')
    return lines


def _f(x, digits=3):
    return 'n/a' if x is None or (isinstance(x, float) and not np.isfinite(x)) else f'{x:.{digits}g}'


def _pct(x):
    return 'n/a' if x is None else f'{100 * x:.1f}%'


def _json_default(o):
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating,)):
        return None if not np.isfinite(o) else float(o)
    if isinstance(o, (np.ndarray,)):
        return o.tolist()
    if isinstance(o, Path):
        return str(o)
    raise TypeError(type(o))


def _clean(o):
    if isinstance(o, float) and not np.isfinite(o):
        return None
    if isinstance(o, dict):
        return {k: _clean(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [_clean(v) for v in o]
    return o


def write(out_dir, result, rows, inputs, protocol_note):
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    result = _clean(json.loads(json.dumps(result, default=_json_default)))
    result['inputs'] = inputs
    result['measurement_protocol'] = protocol_note
    result['reading_bands'] = BANDS
    (out / 'summary.json').write_text(json.dumps(result, indent=2, default=_json_default) + '\n')
    if rows:
        fields = list(dict.fromkeys(k for r in rows for k in r))
        with (out / 'pairs.csv').open('w', newline='') as stream:
            w = csv.DictWriter(stream, fieldnames=fields)
            w.writeheader()
            w.writerows(rows)
    figures = _figures(out, result)
    (out / 'report.html').write_text(_html(result, figures))
    return out


# ----------------------------------------------------------------- figures

def _plt():
    os.environ.setdefault('MPLCONFIGDIR', str(Path.home() / '.cache' / 'matplotlib'))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    return plt


def _figures(out, result):
    plt = _plt()
    heads = result['headline']
    if not heads:
        return []
    paths = []
    if result['study'] == 'convergence':
        x = [h['level'] for h in heads]
        numeric = all(isinstance(v, (int, float)) for v in x)
        xs = x if numeric else list(range(len(x)))
        fig, axes = plt.subplots(1, 3, figsize=(13, 3.8), layout='constrained')
        ax = axes[0]
        ax.plot(xs, [h['path_separation_median_mm'] for h in heads], 'o-', label='typical streamline (median)')
        ax.plot(xs, [h['path_separation_tail_mm'] for h in heads], 's--', label='worst 5% (p95)')
        ax.set(ylabel='separation from reference (mm)', title='Paths: typical vs tail')
        ax = axes[1]
        ax.plot(xs, [h['occupied_jaccard'] for h in heads], 'o-', label='occupied voxels')
        ax.plot(xs, [h['endpoint_jaccard'] for h in heads], 's--', label='endpoints')
        ax.set(ylabel='Jaccard with reference', ylim=(0, 1.02), title='Space')
        ax = axes[2]
        for key, name in [('kappa_w_rel_median_change', 'κw'), ('tau_rel_median_change', 'τ'),
                          ('length_rel_median_change', 'length')]:
            ax.plot(xs, [100 * (h[key] or 0) for h in heads], 'o-', label=name)
        ax.set(ylabel='median shift vs reference (%)', title='MMF geometry')
        for ax in axes:
            ax.set_xlabel(result['parameter'])
            if numeric:
                ax.set_xscale('log')
                ax.invert_xaxis()
            else:
                ax.set_xticks(xs, [str(v) for v in x])
            ax.grid(alpha=.25)
            ax.legend(fontsize=8)
        if numeric:
            axes[0].set_yscale('symlog', linthresh=1e-4)
        fig.suptitle('Convergence toward the finest level (x axis: coarse → fine)')
    else:
        labels = [(h.get('group') + ': ' if h.get('group') else '') + _short(h['a']) + ' vs ' + _short(h['b'])
                  for h in heads]
        y = np.arange(len(heads))
        fig, axes = plt.subplots(1, 3, figsize=(13, 0.9 + 0.45 * len(heads) + 1.5), layout='constrained',
                                 sharey=True)
        tiers = [
            ('Space', [('occupied_jaccard', 'occupied Jaccard'), ('endpoint_jaccard', 'endpoint Jaccard')], (0, 1)),
            ('Per-streamline agreement', [('kappa_w_paired_spearman', 'κw Spearman'),
                                          ('tau_paired_spearman', 'τ Spearman')], (0, 1)),
            ('Distribution agreement', [('kappa_w_rel_median_change', 'κw shift'),
                                        ('tau_rel_median_change', 'τ shift'),
                                        ('length_rel_median_change', 'length shift')], None),
        ]
        for ax, (title, keys, lim) in zip(axes, tiers):
            width = 0.8 / len(keys)
            for k, (key, name) in enumerate(keys):
                vals = [np.nan if h[key] is None else h[key] for h in heads]
                if 'shift' in name:
                    vals = [100 * v for v in vals]
                ax.barh(y + (k - (len(keys) - 1) / 2) * width, vals, width, label=name)
            ax.set_title(title)
            ax.set_xlim(*(lim or (0, None)))
            ax.set_xlabel('%' if lim is None else 'agreement (1 = identical)')
            ax.grid(axis='x', alpha=.25)
            ax.legend(fontsize=8, loc='lower right')
        axes[0].set_yticks(y, labels, fontsize=8)
        axes[0].invert_yaxis()
        fig.suptitle(f"{result['study'].capitalize()}: agreement by tier")
    path = out / 'agreement.png'
    fig.savefig(path, dpi=140)
    plt.close(fig)
    paths.append(path)
    return paths


def _short(label, n=38):
    return label if len(label) <= n else '…' + label[-(n - 1):]


# ----------------------------------------------------------------- html

QUESTION = {
    'convergence': 'Does the output settle as the numerical setting is refined, and does it settle for every streamline or only the typical one?',
    'repeatability': 'How much of the output survives an independent repeat of the same pipeline?',
    'consistency': 'How much does the output depend on the choice of method or setting?',
}


def _html(result, figures):
    heads = result['headline']
    study = result['study']
    parts = [f"<h1>{study.capitalize()} assessment</h1><p class=lede>{QUESTION[study]}</p>"]
    if study == 'convergence':
        ref = result['reference']
        parts.append(f"<p>Knob: <code>{html.escape(str(result['parameter']))}</code>, levels "
                     f"{', '.join(map(str, result['values']))}. Reference: <code>{html.escape(ref['label'])}</code>"
                     + ('' if ref['independent_of_ladder'] else
                        ' (the finest level. This measures <i>self</i>-convergence, not accuracy)') + '.</p>')
        final = heads[-1] if heads else None
        if final:
            parts.append('<h2>Reading (level closest to the reference)</h2><ul>' +
                         ''.join(f'<li>{l}</li>' for l in reading(final)) + '</ul>')
        if result.get('observed_order'):
            orders = ', '.join(f"{o['order']:.2f}" for o in result['observed_order'])
            parts.append(f'<p>Observed order of the median path error between successive levels: {orders}.</p>')
    else:
        for h in heads:
            title = (h.get('group') + ' · ' if h.get('group') else '') + f"{_short(h['a'])} vs {_short(h['b'])}"
            parts.append(f'<h3>{html.escape(title)}</h3><ul>' + ''.join(f'<li>{l}</li>' for l in reading(h)) + '</ul>')
    for f in figures:
        data = base64.b64encode(Path(f).read_bytes()).decode()
        parts.append(f'<figure><img src="data:image/png;base64,{data}" alt="agreement figure"></figure>')
    cols = [('occupied_jaccard', 'Occupied J'), ('endpoint_jaccard', 'Endpoint J'),
            ('path_separation_median_mm', 'Path sep. median (mm)'), ('path_separation_tail_mm', 'Path sep. p95 (mm)'),
            ('kappa_w_rel_median_change', 'κw shift'), ('kappa_w_paired_spearman', 'κw ρ'),
            ('tau_rel_median_change', 'τ shift'), ('tau_paired_spearman', 'τ ρ')]
    first = 'level' if study == 'convergence' else 'pair'
    parts.append('<h2>All comparisons</h2><div class=wrap><table><tr><th>' + first + '</th>' +
                 ''.join(f'<th>{c}</th>' for _, c in cols) + '</tr>')
    for h in heads:
        name = str(h['level']) if study == 'convergence' else f"{_short(h['a'], 28)} / {_short(h['b'], 28)}"
        parts.append('<tr><td>' + html.escape(name) + '</td>' + ''.join(
            f"<td>{_pct(h[k]) if 'shift' in c else _f(h[k])}</td>" for k, c in cols) + '</tr>')
    parts.append('</table></div>')
    parts.append('<h2>How to read this</h2><ul>'
                 '<li><b>Space</b>: voxels the streamlines pass through and end in. This is the only tier that says '
                 'the tractograms are in the same <i>place</i>.</li>'
                 '<li><b>Paths</b>: corresponding streamlines, matched by seed when both runs record seeds (exact), '
                 'otherwise by nearest streamline (heuristic).</li>'
                 '<li><b>MMF geometry</b>: per-streamline 90th percentile of |κ<sub>w</sub>| (curvature, gauge '
                 'invariant, primary), |τ| (torsion) and |w<sub>23</sub>| (frame twist, gauge-dependent diagnostic), '
                 f"in world mm with protocol {html.escape(json.dumps(result['measurement_protocol']))}. "
                 'The shape of a streamline does not tell you whether it is in the right place.</li>'
                 f"<li>The words high/moderate/low are reading conventions ({html.escape(json.dumps(BANDS))}), "
                 'not validated thresholds.</li></ul>')
    parts.append('<p class=small>Full numbers: <code>summary.json</code> (every tier, every comparison) and '
                 '<code>pairs.csv</code> (one row per matched streamline pair).</p>')
    css = ('body{font:15px/1.6 system-ui,sans-serif;max-width:1040px;margin:auto;padding:24px 16px;color:#1d2733;'
           'background:#fff}h1{margin-bottom:0}.lede{font-size:18px;color:#46566a}figure{margin:20px 0}'
           'img{max-width:100%}table{border-collapse:collapse;font-size:13px}td,th{padding:6px 10px;'
           'border-bottom:1px solid #dde3ea;text-align:left}th{background:#eef3f8}.wrap{overflow-x:auto}'
           'code{background:#f1f4f7;padding:1px 4px;border-radius:3px}.small{font-size:13px;color:#5a6878}')
    return ('<!doctype html><html lang=en><head><meta charset=utf-8><meta name=viewport '
            'content="width=device-width,initial-scale=1"><title>Tractography ' + study +
            '</title><style>' + css + '</style></head><body>' + ''.join(parts) + '</body></html>')
