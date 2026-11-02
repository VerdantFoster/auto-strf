"""Observation-specific tracking TLE, cached independently of current catalogs."""
import json
from pathlib import Path
import requests
from .doppler import validate_tle


def scientific_metadata(meta):
    # Do not persist API user/observer details, credentials or host information.
    keys=('id','norad_cat_id','ground_station','start','end','station_lat',
          'station_lng','station_alt','tle0','tle1','tle2','tle_source','station_name','waterfall','frequency')
    return {k:meta[k] for k in keys if k in meta}


def get_correction_tle(observation_id, cache_dir=None):
    """Return the TLE stored on this observation, never the latest catalog TLE.

    Existing observation metadata is the immutable local cache. Geometry remains
    scientific input; operator identity and other unnecessary API fields are omitted.
    """
    oid=int(observation_id)
    if oid<=0:raise ValueError('Observation ID must be positive')
    root=Path(cache_dir) if cache_dir is not None else Path(__file__).resolve().parents[1]/'data/waterfall'
    folder=root/str(oid);cache=folder/'metadata.json'
    if cache.exists():
        meta=json.loads(cache.read_text())
    else:
        response=requests.get(f'https://network.satnogs.org/api/observations/{oid}/',timeout=40)
        response.raise_for_status();meta=scientific_metadata(response.json())
        if meta.get('id')!=oid:raise ValueError('Observation ID mismatch')
        validate_tle(meta)
        folder.mkdir(parents=True,exist_ok=True)
        tmp=folder/'metadata.json.part';tmp.write_text(json.dumps(meta,indent=2));tmp.replace(cache)
    if meta.get('id')!=oid:raise ValueError('Cached observation ID mismatch')
    validate_tle(meta)
    return meta['tle1'],meta['tle2']
