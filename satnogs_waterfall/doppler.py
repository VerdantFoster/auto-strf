"""Undo a declared correction exactly once, on observed frequencies (Hz).
The ideal tracking-TLE model is an explicit assumption, not receiver telemetry.
"""
import json
from functools import lru_cache
import numpy as np
from sgp4.api import Satrec, WGS72
C_M_S = 299792458.0


def correction_state(meta):
    # Standard SatNOGS waterfall is tracked. A correction update-rate setting
    # is NOT an on/off flag: WorkerFreq separately sends tracked RF to rigctld.
    # This is a pipeline input convention, not proof of actual tuner telemetry.
    return 'corrected'



def validate_tle(meta):
    try:
        sat = _validated_tle(meta['tle1'],meta['tle2'])
        if sat.satnum != int(meta['norad_cat_id']):
            raise ValueError('Tracking TLE does not match observation satellite')
        return sat
    except (KeyError,TypeError) as exc:
        raise ValueError('Tracking TLE missing or malformed') from exc


@lru_cache(maxsize=64)
def _validated_tle(l1,l2):
    try:
        for n, line in enumerate((l1, l2), 1):
            if not isinstance(line, str) or len(line) != 69 or not line.startswith(f'{n} '):
                raise ValueError('Expected TLE with 69-character lines')
            if sum(int(c) if c.isdigit() else int(c == '-') for c in line[:68]) % 10 != int(line[-1]):
                raise ValueError('Invalid tracking-TLE checksum')
        if l1[2:7] != l2[2:7]:
            raise ValueError('TLE lines have different satellite IDs')
        sat = Satrec.twoline2rv(l1, l2, WGS72)
        if sat.sgp4(sat.jdsatepoch, sat.jdsatepochF)[0]:
            raise ValueError('Tracking TLE cannot propagate')
        return sat
    except (KeyError, TypeError) as exc:
        raise ValueError('Tracking TLE missing or malformed') from exc


def tracking_shift_hz(times, meta, center_hz):
    times = np.asarray(times, dtype=float)
    if times.ndim != 1 or not len(times) or not np.isfinite(times).all():
        raise ValueError('Undoppler requires finite, nonempty MJD times')
    try:
        lat, lon, alt = [float(meta[k]) for k in ('station_lat','station_lng','station_alt')]
    except (KeyError, TypeError, ValueError) as exc:
        raise ValueError('Receiver coordinates missing or invalid') from exc
    if not np.isfinite([lat,lon,alt,center_hz]).all() or abs(lat)>90 or abs(lon)>180 or not -1000<alt<100000 or center_hz<=0:
        raise ValueError('Receiver coordinates or center frequency outside supported range')
    sat=validate_tle(meta)
    errors,r,v=sat.sgp4_array(np.floor(times)+2400000.5,times%1)
    if np.any(errors) or not np.isfinite(r).all() or not np.isfinite(v).all():
        raise ValueError('SGP4 failed at measurement times')
    jd=times+2400000.5;centuries=(jd-2451545.0)/36525
    theta=np.deg2rad((280.46061837+360.98564736629*(jd-2451545)+.000387933*centuries**2-centuries**3/38710000)%360+lon)
    lat=np.deg2rad(lat);e2=(1/298.257223563)*(2-1/298.257223563)
    n=6378.137/np.sqrt(1-e2*np.sin(lat)**2);h=alt/1000
    site=np.c_[(n+h)*np.cos(lat)*np.cos(theta),(n+h)*np.cos(lat)*np.sin(theta),np.full(len(times),(n*(1-e2)+h)*np.sin(lat))]
    site_v=np.c_[-7.2921159e-5*site[:,1],7.2921159e-5*site[:,0],np.zeros(len(times))]
    rho=r-site;distance=np.linalg.norm(rho,axis=1)
    if np.any(distance<1):raise ValueError('Invalid satellite range')
    radial=np.sum((v-site_v)*rho,axis=1)/distance*1000
    return -float(center_hz)*radial/C_M_S


def remove_doppler_correction(times, offset_hz, meta, center_hz):
    """Tabulation Helper inverse, vectorized for calibrated PNG points.

    Upstream remove_doppler_correction(t, freq) uses freq as a baseband
    OFFSET, not absolute RF: f_center + freq - f_center*range_velocity/c.
    Here SGP4 supplies topocentric range velocity in the same sign convention.
    Attribution and implementation differences are documented in README.md.
    """
    offset=np.asarray(offset_hz,dtype=float)
    shift=tracking_shift_hz(times,meta,center_hz)
    if offset.shape!=shift.shape or not np.isfinite(offset).all():
        raise ValueError('Tabulation Helper requires matching finite time/offset arrays')
    return float(center_hz)+offset+shift


def undo_correction(times, recorded_hz, meta, center_hz, state='metadata', evidence=None):
    t,f=np.asarray(times,dtype=float),np.asarray(recorded_hz,dtype=float)
    if t.ndim!=1 or f.shape!=t.shape or len(t)<2 or not np.isfinite(t).all() or not np.isfinite(f).all() or np.any(np.diff(t)<=0):
        raise ValueError('Undoppler requires at least two chronological finite time/frequency points')
    declared=correction_state(meta)
    if state not in ('metadata','uncorrected','corrected'):raise ValueError('Unknown correction state')
    if state!='metadata' and not (isinstance(evidence,str) and evidence.strip()):raise ValueError('Explicit correction override requires evidence/reason')
    effective=declared if state=='metadata' else state
    if effective=='unknown':raise ValueError('Correction state unknown; provide explicit state and evidence')
    restored=remove_doppler_correction(t,f-float(center_hz),meta,center_hz) if effective=='corrected' else f.copy()
    shift=restored-f
    # Inverse only: f_raw = (f_pixel - f_target) + f_target*(1-vr/c).
    # Never feed restored RF back into this function as PNG frequencies.
    result=restored
    if not np.isfinite(result).all() or np.any(result<=0):raise ValueError('Invalid restored RF frequency')
    return result, {'state':effective,'metadata_state':declared,'model':'ideal_continuous_tracking_tle' if effective=='corrected' else 'identity',
        'evidence':evidence or 'Standard SatNOGS tracked-waterfall convention; ideal reconstruction','correction_applied_once':effective=='corrected','min_added_shift_hz':float(shift.min()),'max_added_shift_hz':float(shift.max()),
        'assumption':'No actual tuner/NCO log available; ideal continuous model only' if effective=='corrected' else 'Metadata declaration, not measured tuning telemetry'}


def diagnostics(times, recorded, restored, meta, target_hz):
    sat=validate_tle(meta)
    epoch=sat.jdsatepoch+sat.jdsatepochF-2400000.5
    age=float(np.min(times)-epoch)
    span=float(np.ptp(recorded))
    return {'correction_tle_epoch_mjd':float(epoch),'tle_age_days_at_first_point':age,
        'residual_span_hz':span,'restored_span_hz':float(np.ptp(restored)),
        'residual_max_abs_hz':float(np.max(np.abs(np.asarray(recorded)-target_hz))),
        'restored_max_abs_hz':float(np.max(np.abs(np.asarray(restored)-target_hz))),
        'shape_note':'Small residual with recent TLE is compatible with successful tracking; shape alone does not validate signal identity or orbit',
        'freshness_note':'Epoch age is diagnostic; no rejection threshold is imposed'}
