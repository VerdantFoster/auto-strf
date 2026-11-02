# SPDX-License-Identifier: AGPL-3.0-or-later
# Interactive event workflow adapted from Fabian P. Schmidt's
# SatNOGS Waterfall Tabulation Helper v0.3 (e5baa1f).
"""Manual tabulation on metadata-calibrated PNGs; no automatic ridge selection."""
import json
from pathlib import Path
import numpy as np
from PIL import Image
from .calibration import read_calibration, sha256


from .doppler import correction_state, undo_correction, diagnostics


def export_points(folder, points, *, state="metadata", evidence=None, context=None, source_method="manual_satnogs_waterfall"):
    folder=Path(folder)
    meta,cal=context if context is not None else (json.loads((folder/'metadata.json').read_text()),read_calibration(folder/'waterfall.png'))
    if not isinstance(meta.get('ground_station'),int) or not 0<meta['ground_station']<=9999: raise ValueError('STRF station ID must be 1..9999')
    points=np.asarray(points,dtype=float)
    if points.ndim!=2 or points.shape[1]!=2 or len(points)<2 or not np.isfinite(points).all():
        raise ValueError('Select at least two finite points')
    x,y=points.T
    if np.any((x<cal.left)|(x>cal.right)|(y<cal.top)|(y>cal.bottom)):
        raise ValueError('Point outside calibrated plot')
    t,hz=cal.pixel_to_world(x,y);order=np.argsort(t);t,hz=t[order],hz[order]
    if np.any(np.diff(t)<=0):raise ValueError('Select at most one frequency per timestamp')
    recorded_hz=hz.copy()
    hz,correction=undo_correction(t,hz,meta,cal.center_hz,state,evidence)
    diagnostic=diagnostics(t,recorded_hz,hz,meta,cal.center_hz)
    print(json.dumps({'observation':meta['id'],**diagnostic}))
    dat=folder/'measurements.dat'
    tmp=folder/'measurements.dat.tmp'
    np.savetxt(tmp,np.c_[t,hz,np.ones(len(t)),np.full(len(t),meta['ground_station'])],fmt=['%.10f','%.6f','%.1f','%d'])
    payload={'export_schema':3,'diagnostics':diagnostic,'coordinate_convention':'imshow_integer_pixel_centers','method':source_method,'observation_id':meta['id'],'points_xy':points[order].tolist(),
             'png_sha256':sha256(folder/'waterfall.png'),'metadata_sha256':sha256(folder/'metadata.json'),
             'dat_sha256':sha256(tmp),'calibration':cal.__dict__,'correction_state':correction['state'],'doppler_correction':correction,'recorded_frequency_hz':recorded_hz.tolist(),'restored_frequency_hz':hz.tolist(),'times_mjd':t.tolist(),
             'point_count':len(t),'frequency_span_hz':float(np.ptp(hz)),
             'time_span_s':float(np.ptp(t)*86400),'sigma_floor_hz':max(cal.native_bin_hz,cal.hz_per_pixel),
             'span_gate_500_hz':False,'validation':'Pixel selection is not orbit validation'}
    (folder/'selection.json.tmp').write_text(json.dumps(payload,indent=2))
    tmp.replace(dat)
    (folder/'selection.json.tmp').replace(folder/'selection.json')
    from .mark_waterfalls import mark_iq_on_waterfall
    mark_iq_on_waterfall(folder/'waterfall.png',dat,folder/'metadata.json',
                         folder/'selection.png',calibration=cal)
    return payload


def dialog(folder, *, state="metadata", evidence=None, context=None):
    import matplotlib.pyplot as plt
    folder=Path(folder)
    meta,cal=context if context is not None else (json.loads((folder/'metadata.json').read_text()),read_calibration(folder/'waterfall.png'))
    markers=[]
    saved=folder/'selection.json'
    if saved.exists():
        previous=json.loads(saved.read_text())
        if previous['png_sha256']!=sha256(folder/'waterfall.png') or previous['metadata_sha256']!=sha256(folder/'metadata.json'):raise ValueError('Saved selection input hashes differ')
        markers=previous['points_xy']
    elif (folder/'proposal.json').exists():
        previous=json.loads((folder/'proposal.json').read_text())
        if previous['png_sha256']!=sha256(folder/'waterfall.png') or previous['metadata_sha256']!=sha256(folder/'metadata.json'):raise ValueError('Proposal input hashes differ')
        markers=previous['points_xy']
        print('Loaded automatic proposals: review before pressing f. They are not yet measurements.',flush=True)
    fig,ax=plt.subplots(figsize=(8,10));ax.imshow(Image.open(folder/'waterfall.png'))
    guide,=ax.plot([],[],color='cyan',linewidth=.9,label='Display guide only')
    line,=ax.plot([],[],'ro',markersize=3,markerfacecolor='none')
    ax.legend(loc='upper left',fontsize=8)
    ax.set_title(f"{meta['id']} | click: add | u: undo | c: clear | f: save\nCorrection: {correction_state(meta)}; use toolbar to zoom")
    def update():
        ordered=sorted(markers,key=lambda p:-p[1])
        line.set_data([p[0] for p in ordered],[p[1] for p in ordered])
        guide.set_data([],[])
        if len(ordered)>=2 and len({p[1] for p in ordered})==len(ordered):
            from .mark_waterfalls import continuous_trace
            trace=continuous_trace(ordered,cal);guide.set_data(trace[:,0],trace[:,1])
        fig.canvas.draw_idle()
    def on_mouse_press(event):
        toolbar=getattr(fig.canvas.manager,'toolbar',None)
        if event.button!=1 or event.inaxes!=ax or (toolbar and toolbar.mode):return
        x,y=event.xdata,event.ydata
        if cal.left<=x<=cal.right and cal.top<=y<=cal.bottom:
            markers.append([x,y]);update()
    def on_key(event):
        if event.key=='u' and markers:markers.pop();update()
        elif event.key=='c':markers.clear();update()
        elif event.key=='f':
            try:
                result=export_points(folder,markers,state=state,evidence=evidence,context=(meta,cal));print(f"Saved {result['point_count']} points to {folder/'measurements.dat'}",flush=True)
            except ValueError as exc: print(f'Not saved: {exc}',flush=True)
    fig.canvas.mpl_connect('button_press_event',on_mouse_press);fig.canvas.mpl_connect('key_press_event',on_key)
    update();plt.show()
