"""Read legacy float32 STRF spectrum records; not a PNG pipeline stage."""
from datetime import datetime,timezone
from pathlib import Path
import re
import numpy as np


def parse_header(header_b):
    if len(header_b)!=256:raise ValueError('Truncated spectrum header')
    try:
        text=header_b.decode('ASCII').strip('\x00')
        pattern=r'^HEADER\nUTC_START\s+([^\n]+)\nFREQ\s+([^\n]+) Hz\nBW\s+([^\n]+) Hz\nLENGTH\s+([^\n]+) s\nNCHAN\s+([^\n]+)\nNSUB\s+([^\n]+)\nEND\n$'
        m=re.fullmatch(pattern,text)
        if m is None:raise ValueError('Unsupported spectrum header (only float32 records supported)')
        h=dict(utc_start=datetime.fromisoformat(m[1]).replace(tzinfo=timezone.utc),freq=float(m[2]),bw=float(m[3]),length=float(m[4]),nchan=int(m[5]),nsub=int(m[6]))
        if not np.isfinite([h['freq'],h['bw'],h['length']]).all() or min(h['freq'],h['bw'],h['length'],h['nchan'],h['nsub'])<=0:raise ValueError('Invalid spectrum dimensions')
        return h
    except (UnicodeError,TypeError,OverflowError) as exc:raise ValueError('Invalid spectrum header') from exc


def read_spectrum(path):
    files=sorted(Path(path).glob('*.bin'))
    if not files:raise ValueError('No STRF spectrum .bin files found')
    records=[]
    for file in files:
        with file.open('rb') as stream:
            while True:
                block=stream.read(256)
                if not block:break
                header=parse_header(block)
                z=np.fromfile(stream,dtype=np.float32,count=header['nchan'])
                if len(z)!=header['nchan'] or not np.isfinite(z).all():raise ValueError('Truncated or nonfinite spectrum samples')
                records.append((header,z))
    if not records:raise ValueError('Empty spectrum files')
    records.sort(key=lambda row:row[0]['utc_start'])
    first=records[0][0]
    for i,(h,z) in enumerate(records):
        if any(h[k]!=first[k] for k in ('freq','bw','nchan')):raise ValueError('Spectrum grid changes between records')
        if i and h['utc_start']<=records[i-1][0]['utc_start']:raise ValueError('Duplicate spectrum timestamp')
    return np.column_stack([z for h,z in records]),[h for h,z in records]
