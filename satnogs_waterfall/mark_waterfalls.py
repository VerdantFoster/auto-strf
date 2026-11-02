"""Render only original image coordinates; never project undopplered RF on a corrected PNG."""
import json
from pathlib import Path
import numpy as np
from PIL import Image
from .calibration import read_calibration, sha256


def continuous_trace(points_xy, cal):
    """Full-height display guide; inserted vertices are never measurements."""
    xy=np.asarray(points_xy,dtype=float)
    if xy.ndim!=2 or xy.shape[1]!=2 or len(xy)<2 or not np.isfinite(xy).all():
        raise ValueError('At least two finite points required for a display guide')
    xy=xy[np.argsort(-xy[:,1])]
    if np.any(np.diff(xy[:,1])>=0):raise ValueError('Duplicate image timestamps')
    if np.any((xy[:,0]<cal.left)|(xy[:,0]>cal.right)|(xy[:,1]<cal.top)|(xy[:,1]>cal.bottom)):
        raise ValueError('Image points outside calibrated ROI')
    # Extend nearest selected frequency to each time boundary. This is a
    # visual convention only, not evidence of a carrier in an empty row.
    trace=xy.copy()
    if xy[0,1]<cal.bottom:trace=np.vstack(([xy[0,0],cal.bottom],trace))
    if xy[-1,1]>cal.top:trace=np.vstack((trace,[xy[-1,0],cal.top]))
    return trace


def mark_iq_on_waterfall(png_path, dat_path, meta_path, out_path, *, calibration=None):
    # Keep the legacy signature; input is waterfall data, not I/Q.
    png_path,dat_path,meta_path,out_path=map(Path,(png_path,dat_path,meta_path,out_path))
    cal=calibration if calibration is not None else read_calibration(png_path)
    selection=json.loads((dat_path.parent/'selection.json').read_text())
    for key,path in [('png_sha256',png_path),('metadata_sha256',meta_path),('dat_sha256',dat_path)]:
        if selection.get(key)!=sha256(path):raise ValueError('Stale source or measurements; re-export selection')
    xy=np.asarray(selection.get('points_xy'),dtype=float)
    if xy.ndim!=2 or xy.shape[1]!=2 or len(xy)<2 or not np.isfinite(xy).all():raise ValueError('No finite selected image points')
    x,y=xy.T
    if np.any((x<cal.left)|(x>cal.right)|(y<cal.top)|(y>cal.bottom)):raise ValueError('Image points outside calibrated ROI')
    if out_path.resolve()==png_path.resolve():raise ValueError('Cannot overwrite source PNG')
    from matplotlib.figure import Figure
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    with Image.open(png_path) as image:
        rgb=np.asarray(image.convert('RGB'))
        width,height=image.size
    # Preserve source geometry; the full-height guide is display-only.
    trace=continuous_trace(xy,cal)
    proposal_path=dat_path.parent/"proposal.json"
    stride=3.
    if proposal_path.exists():
        proposal=json.loads(proposal_path.read_text())
        stride=float(proposal.get("parameters",{}).get("stride_rows",3))
    if not np.isfinite(stride) or stride<=0:raise ValueError("Invalid display stride")
    dpi=100
    fig=Figure(figsize=(width/dpi,height/dpi),dpi=dpi)
    FigureCanvasAgg(fig)
    ax=fig.add_axes([0,0,1,1])
    ax.imshow(rgb,interpolation='nearest')
    # Cyan stays continuous over gaps and boundary extensions. Red overlays
    # adjacent detected samples only, so absence of evidence remains visible.
    ax.plot(trace[:,0],trace[:,1],color='#00ffff',linewidth=.9)
    ordered=xy[np.argsort(-xy[:,1])]
    for a,b in zip(ordered[:-1],ordered[1:]):
        if a[1]-b[1]<=1.5*stride:
            ax.plot([a[0],b[0]],[a[1],b[1]],color='red',linewidth=.9)
    ax.plot(x,y,linestyle='none',marker='o',markersize=3,
            markerfacecolor='none',markeredgecolor='red',markeredgewidth=.65)
    ax.text(cal.left+8,cal.top+8,
            'RED: selected points | CYAN: display guide only\nTime increases upward: START bottom -> END top',
            color='white',fontsize=6,va='top',bbox={'facecolor':'black','alpha':.8,'pad':3})
    ax.annotate('END',trace[-1],xytext=(10,-10),textcoords='offset points',
                color='cyan',fontsize=7,va='top')
    ax.annotate('START',trace[0],xytext=(10,10),textcoords='offset points',
                color='cyan',fontsize=7,va='bottom')
    ax.set_xlim(-.5,width-.5);ax.set_ylim(height-.5,-.5);ax.set_axis_off()
    out_path.parent.mkdir(parents=True,exist_ok=True)
    fig.savefig(out_path,dpi=dpi)
    provenance={'observation_id':selection['observation_id'],
                'renderer':'helper scientific annotation; Matplotlib Agg',
                'source_png_sha256':sha256(png_path),'selection_sha256':sha256(dat_path.parent/'selection.json'),
                'measurements_sha256':sha256(dat_path),'output_sha256':sha256(out_path),
                'point_count':len(xy),'selection_method':selection.get('method'),
                'coordinate_basis':'original points_xy; not uncorrected RF projection',
                'marker':'red open circles; continuous cyan guide with red detected segments',
                'display_trace_xy':trace.tolist(),'display_start_mjd':cal.t_min,'display_end_mjd':cal.t_max,
                'gap_policy':'linear connection between selected points; nearest frequency extended to PNG time boundaries',
                'display_vertices_are_measurements':False,'fit_input_modified':False,
                'time_direction':'bottom=start; top=end','dimensions':[width,height]}
    out_path.with_suffix('.provenance.json').write_text(json.dumps(provenance,indent=2))
    return len(xy)


def run(folder):
    folder=Path(folder)
    return mark_iq_on_waterfall(folder/'waterfall.png',folder/'measurements.dat',folder/'metadata.json',folder/'selection.png')

if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('folder')
    run(p.parse_args().folder)
