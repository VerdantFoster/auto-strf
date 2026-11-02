"""Offline regression checks: python -m satnogs_waterfall.sanity_check."""
import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch
import numpy as np
from .calibration import Calibration,read_calibration
from .doppler import validate_tle,undo_correction,remove_doppler_correction,C_M_S,correction_state
from .helper import export_points
from .suggest import track
from . import mark_waterfalls


class PipelineChecks(unittest.TestCase):
    def test_generic_tle_identity(self):
        # Public upstream fixture, independent of the removed research dataset.
        l1='1 52736U 22057E   23036.86389985  .00017304  00000-0  82793-3 0  9996'
        l2='2 52736  97.5231 153.7414 0006358  83.3132 276.8826 15.19025879 38769'
        meta={'norad_cat_id':52736,'tle1':l1,'tle2':l2}
        self.assertEqual(validate_tle(meta).satnum,52736)
        with self.assertRaises(ValueError):
            validate_tle(dict(meta,norad_cat_id=25544))
        with self.assertRaises(ValueError):
            validate_tle(dict(meta,tle2=l2[:-1]+'0'))
        mismatched=l2[:2]+'25544'+l2[7:68]
        mismatched+=str(sum(int(c) if c.isdigit() else int(c=='-') for c in mismatched)%10)
        with self.assertRaises(ValueError):
            validate_tle(dict(meta,tle2=mismatched))

    def test_axis_and_roundtrip(self):
        c=Calibration(10,210,20,420,437400000,437450000,60000,60001,437425000,100,1)
        t,f=c.pixel_to_world([10,210],[420,20])
        np.testing.assert_equal(t,[60000,60001]);np.testing.assert_equal(f,[437400000,437450000])
        x,y=c.world_to_pixel(t,f);np.testing.assert_equal(x,[10,210]);np.testing.assert_equal(y,[420,20])

    def test_tabulation_helper_sign_and_amplitude(self):
        t=np.array([60000,60000.001,60000.002]);f0=437425000.
        radial=np.array([-7000.,0,7000.]);shift=-f0*radial/C_M_S
        with patch('satnogs_waterfall.doppler.tracking_shift_hz',return_value=shift):
            actual=remove_doppler_correction(t,np.zeros(3),{},f0)
            np.testing.assert_allclose(actual,f0*(1-radial/C_M_S),atol=1e-7,rtol=0)
            self.assertTrue(10000<actual[0]-f0<11000)
            self.assertTrue(-11000<actual[-1]-f0<-10000)
            restored,info=undo_correction(t,np.full(3,f0),{},f0)
            np.testing.assert_equal(restored,actual)

    def test_rate_zero_is_not_raw_flag(self):
        self.assertEqual(correction_state({'client_metadata':json.dumps({'radio':{'parameters':{'doppler-correction-per-sec':'0'}}})}),'corrected')

    def test_explicit_raw_identity(self):
        f=np.array([437420000.,437430000.])
        out,info=undo_correction([60000,60001],f,{},437425000,state='uncorrected',evidence='Synthetic raw control')
        np.testing.assert_equal(out,f)

    def test_bad_time_and_empty(self):
        for t,f in [([],[]),([2,1],[1,1]),([1,1],[1,1]),([1,2],[np.nan,1])]:
            with self.subTest(t=t):
                with self.assertRaises(ValueError):undo_correction(t,f,{},437425000)

    def test_empty_extraction_and_flat(self):
        for roi in [np.empty((0,0)),np.ones((2,2)),np.full((20,20),np.nan)]:
            with self.assertRaises(ValueError):track(roi)
        # No line: explicit empty detection, not a fake ridge.
        self.assertFalse(track(np.zeros((120,100)))[3].any())
        # A straight bright signal must survive; no minimum frequency-span gate.
        roi=np.zeros((120,100));roi[:,49:52]=200
        self.assertGreater(int(track(roi)[3].sum()),20)

    def test_bad_png(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/'bad.png';p.write_bytes(b'corrupt')
            with self.assertRaises(ValueError):read_calibration(p)

    def test_export_order_repeat_and_once(self):
        from PIL import Image
        c=Calibration(0,99,0,99,437420000,437430000,60000,60000.01,437425000,100,1)
        meta={'id':1,'norad_cat_id':25544,'ground_station':1}
        with tempfile.TemporaryDirectory() as tmp:
            folder=Path(tmp);Image.new('RGB',(100,100)).save(folder/'waterfall.png');(folder/'metadata.json').write_text(json.dumps(meta))
            pts=[[50,10],[50,90],[50,50]]
            with patch('satnogs_waterfall.doppler.tracking_shift_hz',return_value=np.array([10000.,0,-10000.])) as model, patch('satnogs_waterfall.helper.diagnostics',return_value={}):
                first=export_points(folder,pts,context=(meta,c));second=export_points(folder,pts,context=(meta,c))
                self.assertEqual(model.call_count,2) # exactly one per complete export
                self.assertEqual(first['dat_sha256'],second['dat_sha256'])
                dat=np.loadtxt(folder/'measurements.dat')
                self.assertTrue(np.all(np.diff(dat[:,0])>0))
                self.assertAlmostEqual(np.ptp(dat[:,1]),20000)
                with self.assertRaises(ValueError):export_points(folder,[],context=(meta,c))
                with self.assertRaises(ValueError):export_points(folder,[[50,10],[51,10]],context=(meta,c))

    def test_native_invalid_inputs(self):
        import subprocess
        binary=Path(__file__).resolve().parents[1]/'rffit-batch'
        if not binary.exists():self.skipTest('Build rffit-batch first')
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)
            for data in ('','bad\nwrong\n','61289.0 nan 1 1\n61289.1 435880000 1 1\n'):
                (p/'in.dat').write_text(data)
                r=subprocess.run([str(binary),str(p/'in.dat'),'missing.tle',str(p/'out.tle'),str(p/'out.csv'),'mamm','100000'],capture_output=True)
                self.assertNotEqual(r.returncode,0);self.assertFalse((p/'out.tle').exists())

    def test_curved_ridge_center_and_signal_gap(self):
        rng=np.random.default_rng(13)
        height,width=300,140
        center=70+22*np.sin(np.linspace(-1.2,1.2,height))
        x=np.arange(width)
        roi=rng.normal(0,.15,(height,width))
        roi+=24*np.exp(-.5*((x[None,:]-center[:,None])/1.2)**2)
        roi[120:160]=0
        rows,columns,strength,keep=track(roi)
        measured=keep & ((rows<110)|(rows>170))
        self.assertGreater(measured.sum(),60)
        self.assertLess(np.median(abs(columns[measured]-center[rows[measured]])),.2)
        self.assertLess(np.max(abs(columns[measured]-center[rows[measured]])),.6)
        self.assertFalse(keep[(rows>128)&(rows<152)].any())

    def test_continuous_guide_preserves_measurements(self):
        from .mark_waterfalls import continuous_trace
        c=Calibration(0,99,0,100,437420000,437430000,60000,60000.01,437425000,100,1)
        points=np.array([[40.,70.],[50.,20.]])
        original=points.copy()
        trace=continuous_trace(points,c)
        np.testing.assert_equal(points,original)
        np.testing.assert_equal(trace,[[40,100],[40,70],[50,20],[50,0]])
        times,_=c.pixel_to_world(trace[:,0],trace[:,1])
        self.assertTrue(np.all(np.diff(times)>0))
        self.assertEqual(times[0],c.t_min);self.assertEqual(times[-1],c.t_max)
        for bad in ([], [[1,2]], [[40,20],[50,20]]):
            with self.assertRaises(ValueError):continuous_trace(bad,c)

    def test_observation_tle_cache(self):
        from .reference import get_correction_tle
        sample={'id':123,'norad_cat_id':25544,'tle1':'first','tle2':'second','observer':'must not persist'}
        with tempfile.TemporaryDirectory() as tmp, patch('satnogs_waterfall.reference.validate_tle'), patch('satnogs_waterfall.reference.requests.get') as request:
            request.return_value.json.return_value=sample
            self.assertEqual(get_correction_tle(123,tmp),('first','second'))
            self.assertEqual(get_correction_tle(123,tmp),('first','second'))
            self.assertEqual(request.call_count,1)
            self.assertNotIn('observer',json.loads((Path(tmp)/'123/metadata.json').read_text()))


if __name__=='__main__':unittest.main()
