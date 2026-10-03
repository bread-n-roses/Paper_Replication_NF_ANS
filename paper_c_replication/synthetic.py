"""Fictional data generated without reading, fitting, or sampling empirical data."""
import json
import numpy as np
import pandas as pd
from .io import fresh_directory

SEED=314159

def generate(path):
    out=fresh_directory(path)
    rng=np.random.default_rng(SEED)
    scores=[];cov=[];bench=[];cube=[];refs=[]
    n=1800
    ids=np.array([f'SYNTH-{i:05d}' for i in range(n)])
    # Three ample invented fields exercise n>=300; a fourth exercises exclusion.
    fields=np.array(['Synthetic field A']*600+['Synthetic field B']*600+['Synthetic field C']*540+['Synthetic small field']*60)
    for year in (2023,2024):
        raw=rng.integers(4,500,n);filtered=raw-rng.binomial(raw,.17)
        raw[-12:]=0;filtered[-12:]=0;filtered[-30:-12]=0
        latent=rng.lognormal(-.25,.85,n)
        benchmark=latent*np.exp(rng.normal(0,.13,n))
        benchmark[::83]=0;benchmark[::97]=np.nan
        base={u:np.ones(n,dtype=bool) for u in ('J','N','OA','L')}
        base['J'][-190:]=False;base['N'][-90:]=False;base['L'][-140:]=False
        in_state={u: b.copy() for u,b in base.items()}
        in_state['J'][1550:1560]=False
        endpoint={}
        for u in base:
            noise=rng.normal(0,.18,n)
            for treatment,counts,shift in [('Raw',raw,.10),('Filtered',filtered,-.03)]:
                ans=np.maximum(0,latent*np.exp(noise+shift+rng.normal(0,.07,n)))
                ans[::79]=0
                mask=in_state[u]
                ans[~mask|(counts==0)]=np.nan
                share=counts/np.sum(counts[mask])
                factor=np.nansum(ans*share)
                ans/=factor
                nf=np.nan_to_num(ans)*share*100
                endpoint[u,treatment]=ans
                for i in np.flatnonzero(mask):
                    scores.append([year,ids[i],u,treatment,nf[i],ans[i],counts[i],bool(i%31==0),bool(raw[i]==0)])
        for i in range(n):
            cov.append([year,ids[i],fields[i],raw[i],filtered[i],int(rng.poisson(raw[i]*1.6)),*[bool(base[u][i]) for u in ('J','N','OA','L')]])
            if base['J'][i]:bench.append([year,ids[i],benchmark[i]])
            if year==2024:
                totals=rng.multinomial(int(raw[i]),np.ones(5)/5)
                zeros=rng.multivariate_hypergeometric(totals,int(raw[i]-filtered[i]))
                refs.extend([[year,ids[i],pub,int(t),int(z)] for pub,t,z in zip(range(2019,2024),totals,zeros)])
        if year==2024:
            for i in np.flatnonzero(in_state['J']):
                a=endpoint['J','Raw'][i]; b=endpoint['J','Filtered'][i]
                for e in (0,1):
                    for t in (0,1):
                        for d in (0,1):
                            val=a if (e,t,d)==(0,0,0) else b if (e,t,d)==(1,1,1) else ((3-e-t-d)*a+(e+t+d)*b)/3+.02*(e-t)**2
                            cube.append([year,ids[i],f'E{e}T{t}A{d}',val])
    from .io import SCHEMAS
    for name,rows in [('scores',scores),('covariates',cov),('benchmark',bench),('cube',cube),('reference_counts',refs)]:
        pd.DataFrame(rows,columns=SCHEMAS[name]).to_csv(out/f'{name}.csv',index=False)
    pd.DataFrame([{'role':role,'classification':label} for role in ('case','control') for label,count in
       [('Crossref-positive',13),('Crossref-zero',9),('indeterminate',5)] for _ in range(count)]).to_csv(out/'crossref_validation.csv',index=False)
    (out/'input_manifest.json').write_text(json.dumps({'mode':'synthetic','seed':SEED,
        'description':'Independent fictional design demonstration. No fitted or resampled empirical values.',
        'not_paper_results':True},indent=2)+'\n',encoding='utf-8')
    return out
