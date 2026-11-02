"""TLE-independent proposals on calibrated PNG; proposals are not reviewed measurements."""
import json
from pathlib import Path
import numpy as np
from PIL import Image
from scipy.ndimage import gaussian_filter,gaussian_filter1d,label
from scipy.spatial import cKDTree
from .calibration import read_calibration,sha256


def track(intensity):
    intensity=np.asarray(intensity,dtype=float)
    if intensity.ndim!=2 or min(intensity.shape)<20 or not np.isfinite(intensity).all():
        raise ValueError("Extraction requires a finite 2D ROI with at least 20 rows and columns")
    smooth=gaussian_filter(intensity,(1.5,.7))
    high=smooth-gaussian_filter1d(smooth,10,axis=1)
    noise=np.maximum(1.,1.4826*np.median(np.abs(high-np.median(high,axis=1)[:,None]),axis=1))
    z=(high/noise[:,None])[::3];n,w=z.shape
    margin=max(8,int(w*.1));emission=np.arcsinh(np.maximum(z,-3));emission[:,:margin]=-100;emission[:,-margin:]=-100
    score=emission[0].copy();back=np.zeros((n,w),dtype=int);cols=np.arange(w)
    for row in range(1,n):
        choices=[]
        for delta in range(-3,4):
            prev=cols+delta;v=score[np.clip(prev,0,w-1)]-.35*delta**2
            v[(prev<0)|(prev>=w)]=-np.inf;choices.append(v)
        choices=np.array(choices);best=np.argmax(choices,axis=0)
        back[row]=np.clip(cols+best-3,0,w-1);score=emission[row]+choices[best,cols]
    path=np.zeros(n,dtype=int);path[-1]=np.argmax(score)
    for row in range(n-1,0,-1):path[row-1]=back[row,path[row]]
    # Refine the ridge to its local intensity maximum, then a parabolic
    # subpixel peak. No TLE or hand-drawn target curve is used.
    refined=path.astype(float)
    strength=np.zeros(n)
    for row,col in enumerate(path):
        lo=max(margin,col-2);hi=min(w-margin,col+3)
        peak=lo+int(np.argmax(z[row,lo:hi]))
        refined[row]=peak;strength[row]=z[row,peak]
        if 0<peak<w-1:
            left,center,right=z[row,peak-1:peak+2]
            curvature=left-2*center+right
            if curvature<0:
                refined[row]+=np.clip(.5*(left-right)/curvature,-.5,.5)
    # Hysteresis: weak continuation is accepted only in a sustained component
    # containing a strong detection. Blank intervals remain unmeasured.
    weak=strength>=2.5
    groups,count=label(weak);keep=np.zeros(n,dtype=bool)
    for group in range(1,count+1):
        mask=groups==group
        if mask.sum()>=8 and strength[mask].max()>=5:keep|=mask
    return np.arange(n)*3,refined,strength,keep



def suggest(folder, *, context=None):
    import matplotlib as mpl
    folder=Path(folder);cal=context[1] if context is not None else read_calibration(folder/'waterfall.png')
    with Image.open(folder/'waterfall.png') as image: rgb=np.asarray(image.convert('RGB'))
    left,right=int(np.ceil(cal.left))+1,int(np.floor(cal.right))-1
    top,bottom=int(np.ceil(cal.top))+1,int(np.floor(cal.bottom))-1
    crop=rgb[top:bottom,left:right]
    if crop.size==0 or min(crop.shape[:2])<20: raise ValueError('Empty or too small waterfall ROI')
    palette=np.asarray(mpl.colormaps['viridis'](np.arange(256)))[:,:3]*255
    distance,index=cKDTree(palette).query(crop.reshape(-1,3))
    if np.mean(distance<8)<.98:raise ValueError('PNG colors do not match expected palette')
    rows,cols,score,keep=track(index.reshape(crop.shape[:2]).astype(float))
    xy=np.c_[left+cols[keep],top+rows[keep]];t,f=cal.pixel_to_world(xy[:,0],xy[:,1])
    order=np.argsort(t);xy,t,f=xy[order],t[order],f[order]
    result={'status':'proposed' if len(xy) else 'no_signal','reason':None if len(xy) else 'No ridge passed contrast and segment-length thresholds','coordinate_convention':'imshow_integer_pixel_centers','method':'automatic_png_proposal_v3_ridge_center','reviewed':False,'points_xy':xy.tolist(),'png_sha256':sha256(folder/'waterfall.png'),'metadata_sha256':sha256(folder/'metadata.json'),
            'n_points':len(xy),'time_span_s':float(np.ptp(t)*86400) if len(t) else 0,'frequency_span_hz':float(np.ptp(f)) if len(f) else 0,
            'parameters':{'stride_rows':3,'weak_contrast':2.5,'strong_contrast':5,'minimum_segment_rows':8,'max_step_pixels':3,'transition_penalty':.35,'emission':'asinh contrast','subpixel':'local peak plus three-point parabolic refinement'},'score_note':'Palette contrast, not physical SNR or identity confidence','tle_used':False,'span_gate_500_hz':False}
    (folder/'proposal.json').write_text(json.dumps(result,indent=2))
    print(json.dumps({k:v for k,v in result.items() if k!='points_xy'},indent=2));return result
