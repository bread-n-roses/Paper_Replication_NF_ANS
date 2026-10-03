"""Verify the schema, integrity and numeric invariants of the ONLY real data shipped."""
from pathlib import Path
import gzip
import hashlib
import io
import json
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]
COLUMNS=['year','journal_id','treatment','nf','ans','article_count','is_dangling','is_isolated']

def audit(root=ROOT):
    folder=root/'data/open';meta=json.loads((folder/'provenance.json').read_text())
    for filename,record,whole in [('leiden_scores.csv.gz',meta,True),
                                 ('norwegian_additional_nodes.csv.gz',meta['additional_n_nodes'],False)]:
        blob=(folder/filename).read_bytes()
        if hashlib.sha256(blob).hexdigest()!=record['sha256']:raise ValueError('Open-data hash mismatch.')
        decoded=gzip.decompress(blob)
        frame=pd.read_csv(io.BytesIO(decoded),dtype={'journal_id':'string'},float_precision='round_trip')
        if list(frame.columns)!=COLUMNS:raise ValueError('Unreviewed data column.')
        if len(frame)!=record['rows']:raise ValueError('Open row inventory mismatch.')
        if frame.duplicated(['year','journal_id','treatment']).any():raise ValueError('Duplicate score row.')
        if frame.journal_id.isna().any():raise ValueError('Missing open identity.')
        if set(frame.year)!={2023,2024} or set(frame.treatment)!={'RAW','REF'}:raise ValueError('Wrong years/treatments.')
        for col in ('nf','ans','article_count'):
            values=frame[col].dropna().to_numpy(float)
            if (values<0).any() or not np.isfinite(values).all():raise ValueError('Invalid open score/count.')
        if frame[['nf','article_count']].isna().any().any():raise ValueError('Missing NF/count.')
        if not frame.ans.notna().equals(frame.article_count.gt(0)):raise ValueError('ANS definition does not match eligible mass.')
        for year,annual in frame.groupby('year'):
            raw=annual.loc[annual.treatment.eq('RAW')].set_index('journal_id').sort_index()
            ref=annual.loc[annual.treatment.eq('REF')].set_index('journal_id').sort_index()
            if not raw.index.equals(ref.index):raise ValueError('RAW/REF state drift.')
            if (ref.article_count>raw.article_count).any():raise ValueError('REF article count exceeds RAW.')
        if whole:
            for _,group in frame.groupby(['year','treatment']):
                if not np.isclose(group.nf.sum(),100,atol=1e-8):raise ValueError('Full L NF mass differs from 100.')
                valid=group.ans.notna()
                expected=group.loc[valid,'nf']/100/(group.loc[valid,'article_count']/group.article_count.sum())
                if not np.allclose(expected,group.loc[valid,'ans'],rtol=1e-10,atol=1e-10):raise ValueError('NF/ANS denominator mismatch.')
    print('Open-data inventory, allowed columns, missingness and score invariants PASS.')
    print('This verifies the declared open files; source lineage review is documented separately.')

if __name__=='__main__':audit()
