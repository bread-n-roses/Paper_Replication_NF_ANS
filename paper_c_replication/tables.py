"""Readable LaTeX counterparts of the paper's numerical tables."""
from pathlib import Path
import numpy as np
import pandas as pd

def latex(frame):
    def value(x):
        if pd.isna(x): return '--'
        if isinstance(x,(float,np.floating)): return f'{x:.3f}'
        escaped={'\\':r'\textbackslash{}','&':r'\&','%':r'\%','$':r'\$','#':r'\#','_':r'\_','{':r'\{','}':r'\}','~':r'\textasciitilde{}','^':r'\textasciicircum{}'}
        return ''.join(escaped.get(c,c) for c in str(x))
    row=lambda xs:' & '.join(value(x) for x in xs)+r' \\'
    return '\n'.join([r'\begin{tabular}{'+'l'*len(frame.columns)+'}',r'\toprule',row(frame.columns),r'\midrule',
                       *[row(xs) for xs in frame.itertuples(index=False,name=None)],r'\bottomrule',r'\end{tabular}',''])

def write(frame,path,mode):
    header='% '+('SYNTHETIC DEMONSTRATION; NOT PAPER RESULTS' if mode=='synthetic' else 'LOCAL EMPIRICAL OUTPUT; DO NOT DISTRIBUTE')+'\n'
    Path(path).write_text(header+latex(frame),encoding='utf-8')

def paired(results,year,scope,universes=('J','N','OA')):
    est=results['contrasts'];ci=results['intervals'];rows=[]
    wanted=[f'Filtered_minus_Raw_{u}' for u in universes]
    wanted += [f'{u}_minus_J_{t}' for u in universes if u!='J' for t in ('Raw','Filtered')]
    for name in wanted:
        point=est.loc[est.year.eq(year)&est.support_scope.eq(scope)&est.contrast_id.eq(name)].iloc[0]
        r={'Contrast':name.replace('Filtered','REF').replace('Raw','RAW').replace('OA','A'),'n':int(point.n)}
        for metric,label in [('mean_delta_ae','Change in MAE'),('mean_abs_movement','Mean absolute ANS movement')]:
            interval=ci.loc[ci.year.eq(year)&ci.support_scope.eq(scope)&ci.contrast_id.eq(name)&ci.statistic.eq(metric)].iloc[0]
            r[label]=point[metric];r[label+' 95% interval']=f'[{interval.ci_lower_95:.3f}, {interval.ci_upper_95:.3f}]'
        rows.append(r)
    return pd.DataFrame(rows)

def paper_tables(r,out,draws,mode):
    table=r['reference_observability'].loc[lambda x:x.year.eq(2024),['publication_year','total','zero','positive','zero_percent']]
    total=table[['total','zero','positive']].sum()
    table=pd.concat([table,pd.DataFrame([{'publication_year':'Total',**total.to_dict(),'zero_percent':100*total['zero']/total['total']}])],ignore_index=True)
    table['positive_percent']=100-table.zero_percent
    write(table,out/'main_table1.tex',mode)
    cov=r['coverage'];current=cov.loc[cov.year.eq(2024)&cov.universe.isin(['J','N','OA'])]
    rows=[]
    for u in ('J','N','OA'):
        a=current.loc[current.universe.eq(u)&current.treatment.eq('Raw')].iloc[0]
        b=current.loc[current.universe.eq(u)&current.treatment.eq('Filtered')].iloc[0]
        rows.append({'Universe':u.replace('OA','A'),'Base set':int(a.base_roster_or_raw_active_n),
            'RAW ANS':int(a.ais_defined_n),'REF ANS':int(b.ais_defined_n),
            'Common sample':int(r['support_sizes'].loc[lambda x:x.year.eq(2024)&x.support_id.eq('s6'),'n'].iloc[0])})
    write(pd.DataFrame(rows),out/'main_table2.tex',mode)
    selection=r['sample_selection'].loc[lambda x:x.year.eq(2024)].set_index('selection_state')
    cols=['n','official_ais_median','raw_articles_median','ref_eligible_articles_median','raw_incoming_eligible_citations_median']
    write(selection[cols].T.rename_axis('Measure').reset_index(),out/'main_table3.tex',mode)
    write(pd.DataFrame([
        ['Citation edges','E','All eligible cited articles/reviews','Targets with observed references'],
        ['Teleportation','T','RAW article shares','REF article shares'],
        ['ANS denominator','D','RAW article shares','REF article shares']],
        columns=['Role','Component','RAW setting','REF setting']),out/'main_table4.tex',mode)
    write(paired(r,2024,'annual_s6'),out/'main_table5.tex',mode)
    cube=r['component_benchmark'].copy();cube['rmse']=np.sqrt(cube.mse)
    cube['cell_id']=cube.cell_id.str.replace('A','D')
    write(cube[['cell_id','n','mae','rmse','spearman']],out/'main_table6_panel_a.tex',mode)
    effect=r['component_effects'].loc[lambda x:x.support_scope.eq('annual_s6')].copy()
    effect['effect']=effect.effect.str.replace('A','D')
    write(effect[['effect','n','mean_absolute_movement']],out/'main_table6_panel_b.tex',mode)
    write(pd.concat([r['common_2024'],r['common_2023']])[['year','universe','treatment','n','mae','rmse','spearman','signed']],out/'main_table7.tex',mode)
    coverage_rows=[]
    for year in (2023,2024):
        for u in ('J','N','OA','L'):
            a=cov.loc[cov.year.eq(year)&cov.universe.eq(u)&cov.treatment.eq('Raw')].iloc[0]
            b=cov.loc[cov.year.eq(year)&cov.universe.eq(u)&cov.treatment.eq('Filtered')].iloc[0]
            coverage_rows.append({'Year':year,'Universe':u.replace('OA','A'),'Base set':int(a.base_roster_or_raw_active_n),
                'RAW network':int(a.frozen_state_n),'RAW ANS':int(a.ais_defined_n),'REF ANS':int(b.ais_defined_n),
                'Loss (%)':100*(a.ais_defined_n-b.ais_defined_n)/a.ais_defined_n,
                'RAW with AIS':int(a.benchmark_evaluable_n),'REF with AIS':int(b.benchmark_evaluable_n)})
    write(pd.DataFrame(coverage_rows),out/'supp_A1_coverage.tex',mode)
    write(r['support_sizes'],out/'supp_A2_supports.tex',mode)
    write(paired(r,2023,'annual_s6'),out/'supp_D1_paired.tex',mode)
    write(paired(r,2024,'maximum_support'),out/'supp_E1_paired.tex',mode)
    write(paired(r,2024,'maximum_support',('J','L')),out/'supp_F1_paired.tex',mode)
    for stem,scope in [('supp_E2_annual','maximum'),('supp_F2_annual','leiden')]:
        write(pd.concat([r[f'{scope}_2024'],r[f'{scope}_2023']])[['year','universe','treatment','n','mae','rmse','spearman','signed']],out/(stem+'.tex'),mode)
