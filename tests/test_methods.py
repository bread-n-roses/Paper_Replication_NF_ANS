import unittest
import numpy as np
import pandas as pd
from paper_c_replication import statistics as s

class Methods(unittest.TestCase):
    def test_paired_error_and_movement_are_distinct(self):
        r=s.named_contrast_statistics([1.,2.,0.],[.5,1.,0.],[1.5,2.,1.])
        self.assertAlmostEqual(r['mean_delta_ae'],0.)
        self.assertAlmostEqual(r['mean_abs_movement'],1.)

    def test_decile_boundary_and_ties(self):
        d=s.official_deciles(np.arange(1,101))
        self.assertEqual(d[89],'D2');self.assertEqual(d[90],'D1')
        tied=s.official_deciles([1,1,1,2,3,4,5,6,7,8])
        self.assertEqual(len(set(tied[:3])),1)

    def test_factorial_effects_known_polynomial(self):
        rows=[]
        for e in (0,1):
            for t in (0,1):
                for d in (0,1):
                    rows.append({'issn_l':'SYNTH-A','cell_id':f'E{e}T{t}A{d}','ais':1+2*e+3*t+5*d+7*e*t})
        _,a=s.component_cube_effects(pd.DataFrame(rows));v=a.set_index('effect').signed_mean
        self.assertAlmostEqual(v['tau_E'],5.5);self.assertAlmostEqual(v['tau_T'],6.5)
        self.assertAlmostEqual(v['tau_A'],5);self.assertAlmostEqual(v['tau_ET'],7)
        self.assertAlmostEqual(v['tau_ETA'],0)

    def test_bootstrap_permutation_invariance_and_constant(self):
        f=pd.DataFrame({'issn_l':['SYNTH-C','SYNTH-A','SYNTH-B','SYNTH-D'],'oa_field':['Y','X','X','Y'],'v':[1.,2.,4.,8.]})
        a=s.PreparedStratifiedBootstrap(f);b=s.PreparedStratifiedBootstrap(f.iloc[::-1])
        x,_=a.mean_interval(a.values('v'),statistic='mean',label='test',draws=79)
        y,_=b.mean_interval(b.values('v'),statistic='mean',label='test',draws=79)
        self.assertEqual(x,y)
        z,_=a.mean_interval(np.ones(4)*3,statistic='mean',label='constant',draws=79)
        self.assertEqual(z['ci_lower_95'],3);self.assertEqual(z['ci_upper_95'],3)

if __name__=='__main__':unittest.main()
