import argparse,json,os,shlex
from pathlib import Path
from datetime import datetime,timezone
import requests
from .calibration import read_calibration
from .helper import dialog,correction_state
ROOT=Path(__file__).resolve().parents[1]
DATA=ROOT/'data/waterfall'


def prepare(folder, *, context=None):
    meta,cal=context if context is not None else (json.loads((folder/'metadata.json').read_text()),read_calibration(folder/'waterfall.png'))
    from .reference import get_correction_tle
    get_correction_tle(meta['id'],cache_dir=folder.parent)
    import numpy as np
    try:
        lat,lon,alt=[float(meta[k]) for k in ('station_lat','station_lng','station_alt')]
        if not np.isfinite([lat,lon,alt]).all() or abs(lat)>90 or abs(lon)>180: raise ValueError('Invalid station geometry')
        if not 0<int(meta['ground_station'])<=9999:raise ValueError('Invalid STRF site ID')
    except (KeyError,TypeError,ValueError) as exc:raise ValueError('Invalid receiver metadata') from exc
    (folder/'seed.tle').write_text(f"NORAD {int(meta['norad_cat_id'])}\n"+meta['tle1']+'\n'+meta['tle2']+'\n')
    (folder/'sites.txt').write_text(f"{int(meta['ground_station']):04d} SN {float(meta['station_lat']):9.5f} {float(meta['station_lng']):9.5f} {int(meta['station_alt']):5d}    SatNOGS {meta['ground_station']}\n")
    report={'observation':meta['id'],'station':meta['ground_station'],'start':meta['start'],'correction':correction_state(meta),'calibration':cal.__dict__,'manual_selection_available':(folder/'selection.json').exists()}
    (folder/'inspection.json').write_text(json.dumps(report,indent=2));return report


def fetch_latest(norad_id, stations):
    DATA.mkdir(parents=True,exist_ok=True)
    for station in stations:
        url='https://network.satnogs.org/api/observations/';params={'norad_cat_id':norad_id,'ground_station':station,'ordering':'-start'};found=None
        for page in range(20):
            if not url:break
            r=requests.get(url,params=params,timeout=40);r.raise_for_status();params=None
            entries=r.json()
            eligible=[m for m in entries if m.get('norad_cat_id')==norad_id and m.get('ground_station')==station and m.get('waterfall') and datetime.fromisoformat(m['end'].replace('Z','+00:00'))<=datetime.now(timezone.utc)]
            if eligible:found=max(eligible,key=lambda m:m['start']);break
            url=r.links.get('next',{}).get('url')
        if not found:raise ValueError(f'No completed waterfall found for station {station} within search limit')
        store_observation(found)


def store_observation(found):
    from .reference import scientific_metadata
    found=scientific_metadata(found)
    from .doppler import validate_tle
    validate_tle(found)
    if not found.get("waterfall"):raise ValueError("Observation has no waterfall")
    DATA.mkdir(parents=True,exist_ok=True)
    folder=DATA/str(found['id']);folder.mkdir(exist_ok=True)
    # Inputs with existing selections are immutable.
    if not (folder/'metadata.json').exists():(folder/'metadata.json').write_text(json.dumps(found,indent=2))
    png=folder/'waterfall.png'
    if not png.exists():
        tmp=folder/'waterfall.part'
        try:
            with requests.get(found['waterfall'],stream=True,timeout=60) as download:
                download.raise_for_status();size=0
                with tmp.open('wb') as f:
                    for chunk in download.iter_content(1024*1024):
                        size+=len(chunk)
                        if size>30*1024*1024:raise ValueError('PNG exceeds 30 MiB')
                        f.write(chunk)
            cal=read_calibration(tmp)
            tmp.rename(png)
        finally:
            if tmp.exists():tmp.unlink()
    if 'cal' not in locals():cal=read_calibration(png)
    print(json.dumps(prepare(folder,context=(json.loads((folder/'metadata.json').read_text()),cal)),indent=2))

def main():
    parser=argparse.ArgumentParser(description='SatNOGS Waterfall Tabulation Helper, metadata-calibrated adaptation')
    sub=parser.add_subparsers(dest='cmd',required=True)
    p=sub.add_parser('fetch-latest')
    p.add_argument('--norad-id',type=int,required=True)
    p.add_argument('--stations',type=int,nargs='+',required=True)
    p=sub.add_parser('fetch');p.add_argument('observation_ids',type=int,nargs='+')
    for name in ('inspect','mark','suggest','export','fit-command'):
        p=sub.add_parser(name);p.add_argument('observation_id',type=int)
        if name in ('mark','export'):
            p.add_argument('--correction',choices=('metadata','uncorrected','corrected'),default='metadata')
            p.add_argument('--evidence',help='Non-private reason/evidence for correction model; do not include secrets')
        if name=='export':p.add_argument('--source',choices=('proposal','selection'),default='selection')
    args=parser.parse_args()
    if args.cmd=='fetch-latest':fetch_latest(args.norad_id,args.stations);return
    if args.cmd=='fetch':
        for oid in args.observation_ids:
            r=requests.get(f'https://network.satnogs.org/api/observations/{oid}/',timeout=40);r.raise_for_status();meta=r.json()
            if meta['id']!=oid:raise ValueError('Observation ID mismatch')
            store_observation(meta)
        return
    folder=DATA/str(args.observation_id)
    context=(json.loads((folder/'metadata.json').read_text()),read_calibration(folder/'waterfall.png'))
    info=prepare(folder,context=context)
    if args.cmd=='inspect':print(json.dumps(info,indent=2))
    elif args.cmd=='suggest':
        from .suggest import suggest
        suggest(folder,context=context)
    elif args.cmd=='mark':
        if not os.environ.get('DISPLAY') and not os.environ.get('WAYLAND_DISPLAY'):raise ValueError('Interactive marking needs a graphical desktop. Run this command in your desktop terminal.')
        dialog(folder,state=args.correction,evidence=args.evidence,context=context)
    elif args.cmd=='export':
        from .helper import export_points
        from .calibration import sha256
        source=json.loads((folder/(args.source+'.json')).read_text())
        if source.get('png_sha256')!=sha256(folder/'waterfall.png') or source.get('metadata_sha256')!=sha256(folder/'metadata.json'):raise ValueError('Proposal/selection input hashes differ')
        export_points(folder,source.get('points_xy'),state=args.correction,evidence=args.evidence,context=context,source_method=('automatic_unreviewed_ridge_center_v3' if source.get('method')=='automatic_png_proposal_v3_ridge_center' else 'automatic_unreviewed') if args.source=='proposal' else source.get('method','imported_selection'))
    else:
        from .calibration import sha256
        selection=json.loads((folder/'selection.json').read_text())
        for key,name in [('png_sha256','waterfall.png'),('metadata_sha256','metadata.json'),('dat_sha256','measurements.dat')]:
            if selection[key]!=sha256(folder/name):raise ValueError('Selection/input checksum mismatch')
        if selection.get('export_schema')!=3:raise ValueError('Legacy export: review correction state and re-export before fitting')
        quote=lambda p:shlex.quote(str(p))
        print(f"cd {quote(folder)} && ST_DATADIR={quote(ROOT)} ST_SITES_TXT=sites.txt ST_COSPAR={info['station']} {quote(ROOT/'rffit')} -d {quote(folder/'measurements.dat')} -i {int(context[0]['norad_cat_id'])} -c {quote(folder/'seed.tle')} -s {info['station']}")

if __name__=='__main__':
    try:main()
    except (ValueError,OSError,KeyError,TypeError,requests.RequestException) as exc:raise SystemExit(str(exc))
