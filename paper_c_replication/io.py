"""Explicit inputs and output isolation; never search the research checkout."""
from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd

PACKAGE = Path(__file__).resolve().parents[1]
SCHEMAS = {
    'scores': ['year','issn_l','universe','cell_id','ef','ais','a_count','is_dangling','is_isolated'],
    'covariates': ['year','issn_l','oa_field','a_raw','a_filtered','raw_incoming_eligible_citations',
                   'in_base_J','in_base_N','in_base_OA','in_base_L'],
    'benchmark': ['year','issn_l','official_ais'],
    'cube': ['year','issn_l','cell_id','ais'],
    'reference_counts': ['year','issn_l','publication_year','total','zero'],
    'crossref_validation': ['role','classification'],
}

def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024*1024),b''): h.update(chunk)
    return h.hexdigest()

def outside_package(path):
    path=Path(path).expanduser().resolve()
    if path==PACKAGE or PACKAGE in path.parents:
        raise ValueError('Working inputs and outputs must be OUTSIDE the distributable package.')
    return path

def fresh_directory(path):
    path=outside_package(path)
    if path.exists() and any(path.iterdir()):
        raise ValueError('Output directory must be new or empty; existing results are never overwritten.')
    path.mkdir(parents=True,exist_ok=True)
    return path

def boolean(s):
    if pd.api.types.is_bool_dtype(s): return s.astype(bool)
    parsed=s.astype('string').str.lower().map({'true':True,'false':False,'1':True,'0':False})
    if parsed.isna().any(): raise ValueError('Boolean column contains a missing or invalid value.')
    return parsed.astype(bool)

def read_inputs(path):
    path=outside_package(path)
    meta=json.loads((path/'input_manifest.json').read_text(encoding='utf-8'))
    if meta.get('mode') not in ('synthetic','empirical'):
        raise ValueError('input_manifest.json must identify synthetic or empirical inputs.')
    frames={}; hashes={}
    for name,columns in SCHEMAS.items():
        file=path/(name+'.csv')
        frame=pd.read_csv(file,dtype={'issn_l':'string'},float_precision='round_trip')
        if set(frame.columns)!=set(columns):
            raise ValueError(f'{name}: columns must exactly match the documented schema.')
        if frame.empty: raise ValueError(f'{name}: empty input.')
        if 'year' in frame and not frame.year.isin([2023,2024]).all():
            raise ValueError(f'{name}: only score years 2023 and 2024 are accepted.')
        if 'issn_l' in frame and (frame.issn_l.isna().any() or frame.issn_l.str.strip().eq('').any()):
            raise ValueError(f'{name}: missing identity.')
        if name!='crossref_validation':
            keys=['year','issn_l']
            if name=='scores': keys+=['universe','cell_id']
            if name=='cube': keys+=['cell_id']
            if name=='reference_counts': keys+=['publication_year']
            if frame.duplicated(keys).any(): raise ValueError(f'{name}: duplicate key.')
        for col in columns:
            if col.startswith('in_base_') or col in ('is_dangling','is_isolated'): frame[col]=boolean(frame[col])
        frames[name]=frame;hashes[file.name]=sha(file)
    if set(frames['scores'].universe)!={'J','N','OA','L'}:
        raise ValueError('Scores must cover J, N, OA and L.')
    if set(frames['scores'].cell_id)!={'Raw','Filtered'}:
        raise ValueError('Endpoint treatments must be Raw and Filtered.')
    if set(frames['cube'].cell_id)!={f'E{e}T{t}A{d}' for e in (0,1) for t in (0,1) for d in (0,1)}:
        raise ValueError('The component input must contain the eight registered cells.')
    if not frames['cube'].year.eq(2024).all(): raise ValueError('Component cells must refer to 2024.')
    if not frames['reference_counts'].year.eq(2024).all() or set(frames['reference_counts'].publication_year)!=set(range(2019,2024)):
        raise ValueError('Reference counts must cover the 2019–2023 cited window for score year 2024.')
    if not frames['crossref_validation'].role.isin(['case','control']).all() or not frames['crossref_validation'].classification.isin(['Crossref-positive','Crossref-zero','indeterminate']).all():
        raise ValueError('Invalid frozen Crossref sample role or classification.')
    if set(frames['crossref_validation'].role)!={'case','control'}: raise ValueError('Both Crossref sample roles are required.')
    for name in ('scores','covariates','benchmark'):
        if set(frames[name].year)!={2023,2024}: raise ValueError(f'{name}: both score years are required.')
    keys=set(zip(frames['covariates'].year,frames['covariates'].issn_l))
    for name in ('scores','benchmark','cube'):
        if not set(zip(frames[name].year,frames[name].issn_l))<=keys:
            raise ValueError(f'{name}: an identity is absent from the covariate union.')
    for name,columns in {'scores':['ef','ais','a_count'],'covariates':['a_raw','a_filtered','raw_incoming_eligible_citations'],
                          'benchmark':['official_ais'],'cube':['ais'],'reference_counts':['total','zero']}.items():
        for col in columns:
            v=pd.to_numeric(frames[name][col],errors='raise')
            if ((v.dropna()<0)|~np.isfinite(v.dropna())).any(): raise ValueError(f'{name}: invalid {col}.')
            if col not in ('ais','official_ais') and v.isna().any(): raise ValueError(f'{name}: missing {col}.')
    if (frames['covariates'].a_filtered>frames['covariates'].a_raw).any(): raise ValueError('REF count exceeds RAW.')
    if (frames['reference_counts'].zero>frames['reference_counts'].total).any(): raise ValueError('Zero-reference count exceeds total.')
    return frames,meta,hashes
