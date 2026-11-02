"""Fit an observation and save diagnostic and full-frame STRF PGPLOT PNGs."""
import argparse
import json
import os
import re
from pathlib import Path
import subprocess
import tempfile
import numpy as np
from PIL import Image
from .calibration import sha256

ROOT = Path(__file__).resolve().parents[1]
OUTPUT_NAMES = ('rffit_candidate.tle', 'rffit_residuals.csv', 'rffit.log',
                'rffit_result.png', 'rffit_full.png', 'rffit_provenance.json')
PARAMETER_INDICATORS = {'fixed': '8', 'ma': '58', 'mamm': '568',
                        'raan': '28', 'maraan': '258', 'mammraan': '2568'}


def export_full_frame(folder, tmp, env, mode, frequency_hz, expected_rms):
    """Render the serialized candidate with the batch RF, without another fit."""
    meta = json.loads((folder / 'metadata.json').read_text())
    result = subprocess.run([
        str(ROOT / 'rffit'), '-d', str(folder / 'measurements.dat'),
        '-c', str(tmp / 'rffit_candidate.tle'), '-i', str(int(meta['norad_cat_id'])),
        '-s', str(int(meta['ground_station'])), '-f', str(frequency_hz),
        '-a', PARAMETER_INDICATORS[mode], '-P', str(tmp / 'rffit_full.png') + '/png'],
        env=env, capture_output=True, text=True, check=False)
    if result.returncode:
        raise ValueError(f'Full-frame export failed (exit {result.returncode}): {result.stderr.strip()}')
    match = re.search(r'^render_rms_hz=(\S+)', result.stdout, re.MULTILINE)
    if not match or not np.isclose(float(match[1]), expected_rms, rtol=0, atol=2e-5):
        raise ValueError('Full-frame residual RMS differs from batch output')
    with Image.open(tmp / 'rffit_full.png') as image:
        if image.format != 'PNG':
            raise ValueError('Full-frame renderer did not produce a PNG')
        image.verify()
    return result.stdout


def run(folder, replace=False, mode="mammraan"):
    if mode not in PARAMETER_INDICATORS:
        raise ValueError(f'Unknown fitting mode: {mode}')
    folder = Path(folder).resolve()
    selection = json.loads((folder / 'selection.json').read_text())
    if selection.get('export_schema') != 3:
        raise ValueError('Re-export the selection using the current undoppler pipeline')
    for key, name in [('dat_sha256', 'measurements.dat'),
                      ('png_sha256', 'waterfall.png'), ('metadata_sha256', 'metadata.json')]:
        if selection.get(key) != sha256(folder / name):
            raise ValueError(f'Stale selection: {name} checksum mismatch')
    if not replace and any((folder / name).exists() for name in OUTPUT_NAMES):
        raise ValueError('Fit outputs already exist; use --replace to regenerate')
    measurements = np.atleast_2d(np.loadtxt(folder / 'measurements.dat'))
    if measurements.shape[1] != 4 or len(measurements) < 3 or not np.isfinite(measurements).all():
        raise ValueError('Fit requires at least three finite four-column measurements')
    if np.any(np.diff(measurements[:, 0]) <= 0):
        raise ValueError('Observation measurements must be strictly chronological')
    input_names = ('measurements.dat', 'seed.tle', 'sites.txt', 'selection.json', 'metadata.json')
    hashes = {name: sha256(folder / name) for name in input_names}
    cutoff = float(measurements[-1, 0] + 1)
    with tempfile.TemporaryDirectory(prefix='.rffit-', dir=folder) as tmp:
        tmp = Path(tmp)
        env = os.environ.copy()
        env['ST_SITES_TXT'] = str(folder / 'sites.txt')
        env['ST_DATADIR'] = str(ROOT)
        result = subprocess.run([
            str(ROOT / 'rffit-batch'), str(folder / 'measurements.dat'),
            str(folder / 'seed.tle'), str(tmp / OUTPUT_NAMES[0]),
            str(tmp / OUTPUT_NAMES[1]), mode, str(cutoff), str(tmp / 'rffit_result.png') + '/png'],
            env=env, capture_output=True, text=True, check=False)
        if result.returncode:
            raise ValueError(f'Native fitting failed (exit {result.returncode}): {result.stderr.strip()}')
        match = re.search(r'^frequency_hz=(\S+)', result.stdout, re.MULTILINE)
        frequency_hz = float(match[1]) if match else float('nan')
        if not np.isfinite(frequency_hz) or frequency_hz <= 0:
            raise ValueError('Native fitting did not report a valid fitted RF frequency')
        a = np.atleast_1d(np.genfromtxt(tmp / OUTPUT_NAMES[1], names=True, delimiter=',', encoding='utf-8', dtype=None))
        if len(a) != len(measurements) or not np.allclose(a['mjd'], measurements[:, 0], rtol=0, atol=1e-10):
            raise ValueError('Native output times differ from measurements')
        e = a['residual_hz']
        if not all(np.isfinite(a[k]).all() for k in ('observed_hz', 'predicted_hz', 'residual_hz', 'elevation_deg')):
            raise ValueError('Nonfinite native residual output')
        if not np.allclose(a['observed_hz'], measurements[:, 1], rtol=0, atol=2e-6):
            raise ValueError('Native frequencies differ from exported measurements')
        if not np.allclose(e, a['observed_hz'] - a['predicted_hz'], rtol=0, atol=2e-6):
            raise ValueError('Inconsistent residual calculation')
        metrics = {'n': len(a), 'rmse_hz': float(np.sqrt(np.mean(e**2))),
                   'mae_hz': float(np.mean(abs(e))), 'p95_absolute_hz': float(np.percentile(abs(e), 95)),
                   'bias_hz': float(e.mean()), 'minimum_elevation_deg': float(a['elevation_deg'].min())}
        render_log = export_full_frame(folder, tmp, env, mode, frequency_hz, metrics['rmse_hz'])
        oid = selection['observation_id']
        reviewed = selection.get('method') in ('reviewed_automatic_png_proposal', 'manual_satnogs_waterfall')
        label = 'Titik telah ditinjau; kandidat eksperimental' if reviewed else 'EKSPERIMENTAL — penandaan belum ditinjau'
        provenance = {'observation_id': oid, 'engine': 'native STRF rffit-batch', 'mode': mode,
                      'free_orbital_elements': {'fixed': [], 'ma':['M0'], 'mamm':['M0','mean_motion'], 'raan':['RAAN'], 'maraan':['M0','RAAN'], 'mammraan':['M0','mean_motion','RAAN']}[mode], 'data_usage': 'all selected points used in fit',
                      'selection_method': selection.get('method'), 'review_label': label,
                      'input_sha256': hashes, 'engine_sha256': sha256(ROOT / 'rffit-batch'),
                      'renderer_sha256': sha256(ROOT / 'rffit'), 'fitted_frequency_hz': frequency_hz,
                      'full_frame': {'file': 'rffit_full.png', 'parameter_indicators': PARAMETER_INDICATORS[mode],
                                     'refitted': False, 'rms_matches_batch': True},
                      'metrics': metrics, 'plot_basis': 'native STRF PGPLOT; serialized candidate TLE predictions at measurement times',
                      'orbit_independently_validated': False}
        (tmp / 'rffit.log').write_text(f'Native STRF; mode={mode}; all selected points used\n' + result.stdout
                                     + '\nFull-frame export (candidate TLE and fitted RF; no refit)\n' + render_log)
        for name in OUTPUT_NAMES[:-1]:
            provenance.setdefault('output_sha256', {})[name] = sha256(tmp / name)
        (tmp / 'rffit_provenance.json').write_text(json.dumps(provenance, indent=2))
        if any(sha256(folder / name) != digest for name, digest in hashes.items()):
            raise ValueError('Input changed during fitting; outputs not published')
        for name in OUTPUT_NAMES:
            (tmp / name).replace(folder / name)
    print(json.dumps({'observation': oid, **metrics}))
    return provenance


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('observation_ids', nargs='+', type=int)
    parser.add_argument('--replace', action='store_true')
    parser.add_argument('--mode', choices=('fixed','ma','mamm','raan','maraan','mammraan'), default='mammraan')
    args = parser.parse_args()
    for oid in args.observation_ids:
        run(ROOT / 'data/waterfall' / str(oid), args.replace, args.mode)


if __name__ == '__main__':
    try:
        main()
    except (ValueError, OSError, KeyError) as exc:
        raise SystemExit(str(exc))
