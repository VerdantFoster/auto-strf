"""Summarize native STRF batch results without changing measurements or fits."""
import json,re,sys
from pathlib import Path
from datetime import datetime,timedelta,timezone
import numpy as np
from sgp4.api import Satrec
from .calibration import sha256


def summarize(folder):
    folder=Path(folder);results={}
    for name in ('early_holdout_ma','early_holdout_mamm','all_ma','all_mamm'):
        lines=(folder/(name+'.tle')).read_text().splitlines()[1:3]
        if len(lines)!=2:raise ValueError('TLE must contain two element lines')
        for line in lines:
            if len(line)!=69 or sum(int(c) if c.isdigit() else int(c=='-') for c in line[:68])%10!=int(line[-1]):raise ValueError('Invalid TLE checksum')
        sat=Satrec.twoline2rv(*lines)
        a=np.genfromtxt(folder/(name+'.csv'),names=True,delimiter=',',dtype=None,encoding='utf-8')
        a=np.atleast_1d(a)
        required=('mjd','site_id','split','observed_hz','predicted_hz','residual_hz','elevation_deg')
        if not a.size or not a.dtype.names or not set(required)<=set(a.dtype.names):raise ValueError('Empty or malformed residual CSV')
        if not all(np.isfinite(a[k]).all() for k in ('mjd','observed_hz','predicted_hz','residual_hz','elevation_deg')):raise ValueError('Nonfinite residual CSV')
        splits={}
        for split in ('training','holdout'):
            s=a[a['split']==split]
            if len(s)==0:continue
            e=s['residual_hz']
            if not np.allclose(e,s['observed_hz']-s['predicted_hz'],atol=2e-6,rtol=0):raise ValueError('Residual differs from observed minus predicted')
            for t in s['mjd']:
                if sat.sgp4(2400000.5+np.floor(t),t%1)[0]:raise ValueError('Candidate SGP4 failed')
            splits[split]={'n':len(s),'rmse_hz':float(np.sqrt(np.mean(e**2))),'mae_hz':float(np.mean(np.abs(e))),'bias_hz':float(np.mean(e)),'std_hz':float(np.std(e)),'p95_absolute_hz':float(np.percentile(np.abs(e),95)),'min_elevation_deg':float(np.min(s['elevation_deg']))}
        match=re.search(r'frequency_hz=([\d.eE+-]+)',(folder/(name+'.log')).read_text())
        if not match:raise ValueError('Missing fitted frequency in log')
        freq=float(match.group(1))
        if not np.isfinite(freq) or freq<=0:raise ValueError('Invalid fitted frequency')
        results[name]={'metrics':splits,'fitted_frequency_hz':freq,'tle_sha256':sha256(folder/(name+'.tle')),'physical_check':'below_horizon' if np.min(a['elevation_deg'])<0 else 'above_horizon_at_selected_points'}
    seed=(folder/'seed.tle').read_text().splitlines()[1];yy=int(seed[18:20]);year=(2000 if yy<57 else 1900)+yy;epoch=datetime(year,1,1,tzinfo=timezone.utc)+timedelta(days=float(seed[20:32])-1)
    report={'status':'experimental_fit_completed','evaluation_basis':'Measured waterfall frequencies; no catalog TLE comparison','reference_tle_validation':False,'seed_role':'Initialization only; epoch preserved','physical_diagnostics_are_catalog_comparison':False,'seed_epoch_utc':epoch.isoformat(),'review':'Review provenance belongs to each selection; this summary does not grant approval','assumption':'Use selection.json correction provenance; actual receiver tuning not independently verified','evaluation_role':'Development temporal split; already inspected observations, not new blind validation','split_basis':'training/holdout labels from supplied native residual CSV; no hardcoded observation IDs','native_frequency_model':'One common fitted RF frequency; not separate per-observation offsets','variants':results,
            'input_sha256':{name:sha256(folder/name) for name in ['combined.dat','seed.tle','sites.txt']}}
    (folder/'report.json').write_text(json.dumps(report,indent=2))
    print(json.dumps(report,indent=2))

if __name__=='__main__':summarize(sys.argv[1])
