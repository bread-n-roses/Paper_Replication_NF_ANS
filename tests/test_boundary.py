import hashlib
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
import pandas as pd
import numpy as np
from paper_c_replication.io import PACKAGE,outside_package,read_inputs
from paper_c_replication.synthetic import generate

spec=importlib.util.spec_from_file_location('release_builder',PACKAGE/'tools/build_release.py')
release=importlib.util.module_from_spec(spec);spec.loader.exec_module(release)

class Boundary(unittest.TestCase):
    def test_outputs_cannot_be_inside_package(self):
        with self.assertRaises(ValueError):outside_package(PACKAGE/'data/private')

    def test_archive_requires_reviewed_bytes_and_no_extra_files(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);p=root/'readme.txt';p.write_text('safe')
            record={'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'bytes':len(p.read_bytes())}
            (root/'release_inventory.json').write_text(json.dumps({'files':{'readme.txt':record}}))
            release.check(root)
            p.write_text('changed')
            with self.assertRaises(ValueError):release.check(root)
            p.write_text('safe');(root/'private.csv').write_text('SYNTHETIC CANARY')
            with self.assertRaises(ValueError):release.check(root)

    def test_generator_does_not_read_empirical_inputs(self):
        with tempfile.TemporaryDirectory() as d:
            with patch.object(pd,'read_csv',side_effect=AssertionError('Generator attempted a data read')),patch.object(pd,'read_parquet',side_effect=AssertionError('Generator attempted a data read')):
                out=generate(Path(d)/'inputs')
            f,meta,_=read_inputs(out)
            self.assertEqual(meta['mode'],'synthetic')
            for name,frame in f.items():
                if 'issn_l' in frame:self.assertTrue(frame.issn_l.str.startswith('SYNTH-').all())
            self.assertEqual(len(f['cube'].cell_id.unique()),8)
            self.assertTrue(f['scores'].ais.isna().any())
            self.assertTrue(f['scores'].ais.eq(0).any())

    def test_duplicate_benchmark_key_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            out=generate(Path(d)/'inputs');f=pd.read_csv(out/'benchmark.csv')
            pd.concat([f,f.iloc[:1]]).to_csv(out/'benchmark.csv',index=False)
            with self.assertRaises(ValueError):read_inputs(out)

    def test_csv_round_trip_preserves_near_ties(self):
        with tempfile.TemporaryDirectory() as d:
            out=generate(Path(d)/'inputs');f=pd.read_csv(out/'benchmark.csv')
            values=np.array([.7,np.nextafter(.7,1.),np.nextafter(.7,0.)])
            f.loc[:2,'official_ais']=values;f.to_csv(out/'benchmark.csv',index=False)
            frames,_,_=read_inputs(out)
            np.testing.assert_array_equal(frames['benchmark'].official_ais.iloc[:3],values)

if __name__=='__main__':unittest.main()
