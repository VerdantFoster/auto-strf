"""Plot legacy STRF records using actual timestamps and full bandwidth."""
from datetime import timedelta
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates


def plot_spectrum(z,headers):
    z=np.asarray(z)
    if z.ndim!=2 or not z.size or z.shape[1]!=len(headers) or not np.isfinite(z).all():raise ValueError('Invalid spectrum matrix/headers')
    fig,ax=plt.subplots(figsize=(12,6))
    mean,std=float(z.mean()),float(z.std())
    # Individual time cells preserve gaps instead of stretching every record to 1 s.
    bw=headers[0]['bw'];edges=np.linspace(-bw/2,bw/2,z.shape[0]+1)*1e-3
    for i,h in enumerate(headers):
        if h['length']<=0 or h['bw']!=bw or (i and h['utc_start']<=headers[i-1]['utc_start']):raise ValueError('Invalid spectrum time/frequency grid')
        t=[mdates.date2num(h['utc_start']),mdates.date2num(h['utc_start']+timedelta(seconds=h['length']))]
        ax.pcolormesh(t,edges,z[:,i:i+1],vmin=mean-2*std,vmax=mean+6*std,cmap='viridis',shading='flat')
    ax.xaxis_date();ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M:%S'))
    ax.set(xlabel='Time (UTC)',ylabel=f"Frequency offset (kHz); center {headers[0]['freq']*1e-6:g} MHz")
    fig.autofmt_xdate();fig.tight_layout();plt.show()
