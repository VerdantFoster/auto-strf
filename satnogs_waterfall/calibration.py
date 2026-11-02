"""Calibrated SatNOGS PNG coordinates. Never infer bandwidth from SDR rate."""
from dataclasses import dataclass
import hashlib
import json
from pathlib import Path

from matplotlib.figure import Figure
from matplotlib.gridspec import GridSpec
import numpy as np
from PIL import Image

from datetime import datetime, timezone

def iso_to_mjd(value):
    dt = datetime.fromisoformat(value.replace("Z", "+00:00"))
    if dt.tzinfo is None: raise ValueError("Timestamp must include timezone")
    return (dt - datetime(1858, 11, 17, tzinfo=timezone.utc)).total_seconds()/86400



def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


@dataclass
class Calibration:
    left: float
    right: float
    top: float
    bottom: float
    f_min: float
    f_max: float
    t_min: float
    t_max: float
    center_hz: float
    native_bin_hz: float
    native_row_s: float

    def pixel_to_world(self, x, y):
        return (self.t_min + (self.bottom - np.asarray(y)) / (self.bottom - self.top)
                * (self.t_max - self.t_min),
                self.f_min + (np.asarray(x) - self.left) / (self.right - self.left)
                * (self.f_max - self.f_min))

    def world_to_pixel(self, mjd, hz):
        return (self.left + (np.asarray(hz) - self.f_min) / (self.f_max - self.f_min)
                * (self.right - self.left),
                self.bottom - (np.asarray(mjd) - self.t_min) / (self.t_max - self.t_min)
                * (self.bottom - self.top))

    @property
    def hz_per_pixel(self):
        return (self.f_max - self.f_min) / (self.right - self.left)

    @property
    def seconds_per_pixel(self):
        return (self.t_max - self.t_min) * 86400 / (self.bottom - self.top)


def _read_calibration(path):
    with Image.open(path) as img:
        img.load()
        if img.format != "PNG": raise ValueError("Expected PNG input")
        width, height = img.size
        try:
            plot = json.loads(img.info['satnogs:wf-plot'])
            data = json.loads(img.info['satnogs:wf-dat'])
        except (KeyError, ValueError) as exc:
            raise ValueError('PNG lacks SatNOGS calibration metadata; cannot measure RF/time') from exc
    fig = Figure(figsize=plot['figsize'])
    # SatNOGS renderer uses Matplotlib's default subplot spacing, 0.2.
    spec = dict(plot['gridspec'])
    spec.setdefault('wspace', .2)
    spec.setdefault('hspace', .2)
    box = GridSpec(**spec)[0].get_position(fig)
    center = float(data['center_freq'])
    t_ref = iso_to_mjd(data['timestamp'])
    t0, t1 = [t_ref + float(t) / 86400 for t in plot['ylim_s']]
    f0, f1 = [center + float(f) * 1000 for f in plot['xlim_kHz']]
    sr, nchan = float(data['samp_rate']), int(data['nchan'])
    if sr <= 0 or nchan <= 0 or float(data['nfft_per_row']) <= 0:
        raise ValueError('Sample rate, channel count, and integrations must be positive')
    cal = Calibration(box.x0 * width-.5, box.x1 * width-.5, (1-box.y1)*height-.5,
                      (1-box.y0)*height-.5, f0, f1, t0, t1, center,
                      sr/nchan, float(data['nfft_per_row'])*nchan/sr)
    vals = np.array(list(cal.__dict__.values()))
    if not np.all(np.isfinite(vals)) or not (0 <= cal.left < cal.right <= width
            and 0 <= cal.top < cal.bottom <= height and f0 < f1 and t0 < t1
            and sr > 0 and nchan > 0 and center > 0):
        raise ValueError('Invalid PNG calibration')
    return cal



def read_calibration(path):
    try:
        return _read_calibration(path)
    except (OSError, KeyError, TypeError, ValueError, ZeroDivisionError, OverflowError) as exc:
        raise ValueError(f'Invalid or unsupported waterfall PNG calibration: {exc}') from exc
