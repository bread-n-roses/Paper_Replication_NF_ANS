"""All downstream tables, figures and quantitative summaries from frozen scores."""
from pathlib import Path
import json
import importlib.metadata
import numpy as np
import pandas as pd
from . import statistics as st, contrasts as ct, plots, supplement as supp
from .io import fresh_directory, read_inputs

UNIVERSES=('J','N','OA','L')
TREATMENTS=('Raw','Filtered')

def panel_for(year,frames):
    panel=frames['covariates'].loc[frames['covariates'].year.eq(year)].copy()
    benchmark=frames['benchmark'].loc[frames['benchmark'].year.eq(year),['issn_l','official_ais']]
    panel=panel.merge(benchmark,on='issn_l',how='left',validate='one_to_one')
    panel['in_jcr_one_to_one']=panel.issn_l.isin(benchmark.issn_l)
    panel['official_ais_defined']=panel.official_ais.notna()
    scores=frames['scores'].loc[frames['scores'].year.eq(year)]
    for u in UNIVERSES:
        raw_ids=set(scores.loc[scores.universe.eq(u)&scores.cell_id.eq('Raw'),'issn_l'])
        ref_ids=set(scores.loc[scores.universe.eq(u)&scores.cell_id.eq('Filtered'),'issn_l'])
        if not raw_ids or raw_ids!=ref_ids: raise ValueError(f'{year} {u}: RAW and REF states must be identical and nonempty.')
        if not raw_ids<=set(panel.loc[panel[f'in_base_{u}'],'issn_l']): raise ValueError('Scoring state lies outside base universe.')
        for t in TREATMENTS:
            cell=f'{t}_{u}';short=f'{u.lower()}__{t.lower()}'
            sub=scores.loc[scores.universe.eq(u)&scores.cell_id.eq(t),['issn_l','ais','ef']].rename(columns={'ais':ct._score_key('ais',cell),'ef':ct._score_key('ef',cell)})
            panel=panel.merge(sub,on='issn_l',how='left',validate='one_to_one')
            panel[f'in_state_{short}']=panel[ct._score_key('ef',cell)].notna()
            panel[f'ais_defined_{short}']=panel[ct._score_key('ais',cell)].notna()
    required=[ct._score_key('ais',c) for c in st.PRIMARY_CELLS]
    panel['in_support_s6']=panel[['official_ais',*required]].notna().all(axis=1)
    for spec in ct.CONTRAST_SPECS:
        panel['in_support_'+spec.maximum_support_id]=panel[['official_ais',*[ct._score_key('ais',c) for c in spec.cells]]].notna().all(axis=1)
    return panel.sort_values('issn_l').reset_index(drop=True)

def coverage(panels,frames):
    rows=[];sizes=[]
    for year,panel in panels.items():
        for c in panel.columns:
            if c.startswith('in_support_'):sizes.append({'year':year,'support_id':c.removeprefix('in_support_'),'n':int(panel[c].sum())})
        for u in UNIVERSES:
            for t in TREATMENTS:
                state=panel[f'in_state_{u.lower()}__{t.lower()}'];defined=panel[f'ais_defined_{u.lower()}__{t.lower()}']
                cell=frames['scores'].loc[frames['scores'].year.eq(year)&frames['scores'].universe.eq(u)&frames['scores'].cell_id.eq(t)]
                base=int(panel[f'in_base_{u}'].sum())
                rows.append({'year':year,'universe':u,'treatment':t,'base_roster_or_raw_active_n':base,
                    'frozen_state_n':int(state.sum()),'ais_defined_n':int(defined.sum()),
                    'ais_undefined_n':int((state&~defined).sum()),'base_outside_state_n':base-int(state.sum()),
                    'ais_zero_n':int(panel[ct._score_key('ais',f'{t}_{u}')].eq(0).sum()),
                    'benchmark_evaluable_n':int((defined&panel.official_ais_defined).sum()),
                    'dangling_nodes':int(cell.is_dangling.sum()),'isolated_nodes':int(cell.is_isolated.sum())})
    return pd.DataFrame(rows),pd.DataFrame(sizes)

def contrasts(panels,draws):
    estimates=[];intervals=[]
    for year,panel in panels.items():
        prepared={}
        for spec in ct.CONTRAST_SPECS:
            for scope,support in ct._contrast_supports(spec):
                sample=panel.loc[ct._support_mask(panel,support)]
                values=ct._contrast_statistic(spec,sample)
                meta={'year':year,'contrast_id':spec.contrast_id,'support_scope':scope,'support_id':support,'n':len(sample)}
                estimates.append({**meta,**values})
                # These are exactly the 48 scalar intervals printed in the accepted paper/supplement.
                if spec.kind!='pair' or (year==2023 and scope!='annual_s6'):continue
                if support not in prepared:prepared[support]=st.PreparedStratifiedBootstrap(sample)
                p=prepared[support];a,b=[p.values(ct._score_key('ais',c)) for c in spec.cells];c=p.values('official_ais')
                for name,vector in [('mean_delta_ae',np.abs(b-c)-np.abs(a-c)),('mean_abs_movement',np.abs(b-a))]:
                    label=f'{year}|{support}|{spec.contrast_id}|{name}'
                    record,_=p.mean_interval(vector,statistic=name,label=label,draws=draws)
                    intervals.append({**meta,**record})
        print(f'Completed {year} contrasts',flush=True)
    return pd.DataFrame(estimates),pd.DataFrame(intervals)

def descriptors(panels,frames):
    selection=[];reference=[];gaps=[];field_inventory=[]
    for year,panel in panels.items():
        eligible=panel.loc[panel.official_ais_defined&panel.in_jcr_one_to_one]
        for flag,label in [(True,'S6'),(False,'not_S6')]:
            group=eligible.loc[eligible.in_support_s6.eq(flag)]
            r={'year':year,'selection_state':label,'n':len(group)}
            for col,prefix in [('official_ais','official_ais'),('a_raw','raw_articles'),('a_filtered','ref_eligible_articles'),('raw_incoming_eligible_citations','raw_incoming_eligible_citations')]:
                for q,key in [(.5,'median'),(.25,'q25'),(.75,'q75')]:r[f'{prefix}_{key}']=group[col].quantile(q)
            selection.append(r)
        common=panel.loc[panel.in_support_s6].copy();common['decile']=st.official_deciles(common.official_ais)
        for kind,groups in [('universe',[(u,panel.loc[panel[f'in_state_{u.lower()}__raw']]) for u in UNIVERSES]),
                            ('decile',list(common.groupby('decile',sort=False)))]:
            for label,group in groups:
                total=float(group.a_raw.sum());zero=float((group.a_raw-group.a_filtered).sum())
                gaps.append({'year':year,'kind':kind,'group':label,'total':total,'zero':zero,'zero_percent':100*zero/total if total else np.nan})
        for field,n in common.oa_field.fillna('Unknown').value_counts().items():
            field_inventory.append({'year':year,'oa_field':field,'n':n,'displayed':n>=300})
        counts=frames['reference_counts'];counts=counts.loc[counts.year.eq(year)&counts.issn_l.isin(common.issn_l)]
        if not counts.empty:
            annual=counts.groupby('publication_year')[['total','zero']].sum().reset_index()
            annual['positive']=annual.total-annual.zero;annual['zero_percent']=100*annual.zero/annual.total;annual['year']=year
            reference.append(annual)
    cross=frames['crossref_validation'].groupby(['role','classification']).size().rename('n').reset_index()
    cross['within_role_percent']=100*cross.n/cross.groupby('role').n.transform('sum')
    return {'sample_selection':pd.DataFrame(selection),'reference_observability':pd.concat(reference,ignore_index=True),
            'reference_gaps':pd.DataFrame(gaps),'field_inventory':pd.DataFrame(field_inventory),'crossref_validation':cross}

def save_frame(out,name,frame):
    frame.to_csv(out/'tables'/f'{name}.csv',index=False)
    from .tables import latex
    (out/'tables'/f'{name}.tex').write_text(latex(frame),encoding='utf-8')

def run(inputs,output,draws=4999,figures=True):
    frames,meta,hashes=read_inputs(inputs)
    if meta['mode']=='empirical' and draws!=4999:raise ValueError('Empirical runs require the registered 4,999 draws.')
    out=fresh_directory(output);(out/'tables').mkdir();(out/'figures').mkdir()
    (out/'DO_NOT_DISTRIBUTE.txt').write_text('Local analysis output. Empirical results may reveal licensed values or membership. Never add this directory to the public package.\n')
    panels={y:panel_for(y,frames) for y in (2023,2024)}
    cov,sizes=coverage(panels,frames)
    results={'coverage':cov,'support_sizes':sizes,**descriptors(panels,frames)}
    common={y:p.loc[p.in_support_s6].copy() for y,p in panels.items()}
    decile,field=plots.prepare_group_error_tables(common)
    results.update(decile_profiles=decile,field_profiles=field)
    estimates,intervals=contrasts(panels,draws)
    results.update(contrasts=estimates,intervals=intervals)
    effects,cube_stats=ct._cube_outputs(frames['cube'],panels[2024])
    results.update(component_effects=effects,component_benchmark=cube_stats)
    for year,panel in panels.items():
        for scope,universes,supports in [
            ('common',('J','N','A'),{'J':'s6','N':'s6','A':'s6'}),
            ('maximum',('J','N','A'),{'J':'max_treatment_J','N':'max_treatment_N','A':'max_treatment_OA'}),
            ('leiden',('J','L'),{'J':'max_rectangle_J_L','L':'max_rectangle_J_L'})]:
            specs=supp.specs(universes,supports);summary=supp.summary(panel,specs);summary['year']=year
            results[f'{scope}_{year}']=summary
            for kind in ('decile','field'):results[f'{scope}_{kind}_{year}']=supp.profile(panel,specs,kind=='field')
            if figures and ((scope=='common' and year==2023) or (scope!='common' and year==2024)):
                prefix='replication' if scope=='common' and year==2023 else scope
                supp.headline(summary,out/'figures'/f'{prefix}_headline_{year}.png',year)
                if scope=='common':
                    plots.plot_six_panel_scatter(common[year],out/'figures'/f'ais_c_vs_ais_o_six_panel_{year}.png',year=year)
                else:
                    supp.scatter(panel,specs,out/'figures'/f'{prefix}_scatter_{year}.png',year)
                for kind in ('decile','field'):supp.heatmap(results[f'{scope}_{kind}_{year}'],specs,out/'figures'/f'{prefix}_{kind}_{year}.png',year,kind)
    if figures:
        headline=results['common_2024'].rename(columns={'signed':'signed_mean_error'})
        plots.plot_headline_connected(headline,out/'figures/headline_metrics_connected_2024.png',year=2024)
        plots.plot_six_panel_scatter(common[2024],out/'figures/ais_c_vs_ais_o_six_panel_2024.png',year=2024)
        for metric,prefix in [('scaled_mae','scaled_mae'),('nrmse','nrmse')]:
            ceiling=plots.shared_heatmap_ceiling(decile,field,metric=metric,year=2024)
            plots.plot_decile_heatmap(decile,out/'figures'/f'{prefix}_decile_heatmap_2024.png',metric=metric,year=2024,color_scale_max=ceiling)
            plots.plot_field_heatmap(field,out/'figures'/f'{prefix}_field_heatmap_2024.png',metric=metric,year=2024,color_scale_max=ceiling,expected_fields=None)
    for name,frame in results.items():save_frame(out,name,frame)
    from .tables import paper_tables
    paper_tables(results,out/'tables',draws,meta['mode'])
    if meta['mode']=='synthetic' and figures:
        from PIL import Image,ImageDraw,ImageFont
        for path in (out/'figures').glob('*.png'):
            img=Image.open(path).convert('RGB')
            labelled=Image.new('RGB',(img.width,img.height+70),'white');labelled.paste(img,(0,70))
            draw=ImageDraw.Draw(labelled)
            draw.text((20,20),'SYNTHETIC DEMONSTRATION - NOT PAPER RESULTS',fill='red',font=ImageFont.load_default(size=24))
            labelled.save(path)
    report={'mode':meta['mode'],'bootstrap_draws':draws,'bootstrap_base_seed':st.BASE_SEED,
        'input_sha256':hashes,'packages':{n:importlib.metadata.version(n) for n in ('numpy','pandas','matplotlib','scipy')},
        'table_count':len(results),'figures':len(list((out/'figures').glob('*.png'))),
        'empirical_results_validated':False,'notice':'Running successfully does not itself establish agreement with the accepted manuscript.'}
    (out/'run_report.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(f'Completed {meta["mode"]} analysis. Outputs are local; do not distribute empirical outputs.',flush=True)
    return results,report
