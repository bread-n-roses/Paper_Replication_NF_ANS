# Result-to-code map

Outputs below are created locally under `outputs/tables` and `outputs/figures`.
No empirical version of these outputs is included in the public archive.

| Paper result | Output | Computation |
|---|---|---|
| Main Table 1, reference observability | `main_table1.tex`, `reference_observability.csv` | Sum journal/publication-year counts on derived common support |
| Main Table 2, coverage | `main_table2.tex`, `coverage.csv` | Base/state membership and defined-score counts |
| Main Table 3, sample composition | `main_table3.tex`, `sample_selection.csv` | Within/outside common support among matched official-AIS records |
| Main Table 4, component design | `main_table4.tex` | Three-channel design; no empirical input |
| Main Table 5, paired comparisons | `main_table5.tex`, `contrasts.csv`, `intervals.csv` | Paired point estimates and field-stratified intervals |
| Main Table 6, component experiment | `main_table6_panel_a.tex`, `main_table6_panel_b.tex`, `component_*.csv` | Eight-cell benchmark statistics and seven factorial effects |
| Main Table 7, annual comparison | `main_table7.tex`, `common_2023.csv`, `common_2024.csv` | Same six-cell statistics in each year |
| Main Figure 1 | `headline_metrics_connected_2024.png` | MAE, RMSE, Spearman and signed error |
| Main Figure 2 | `ais_c_vs_ais_o_six_panel_2024.png` | Matched scores and full-range insets; filename retained from manuscript |
| Main Figure 3 | `scaled_mae_decile_heatmap_2024.png` | Group-scaled MAE by AIS decile |
| Main Figure 4 | `scaled_mae_field_heatmap_2024.png` | Group-scaled MAE by rolling field |
| Supplement A.1–A.2 | `supp_A1_coverage.tex`, `supp_A2_supports.tex` | Coverage and comparison-support accounting |
| Supplement B | Upstream documentation / bundled full L scores / Opindx | NF/ANS construction, identity decisions and solver checks; not recomputed |
| Supplement C.1–C.2 figures | `nrmse_decile_heatmap_2024.png`, `nrmse_field_heatmap_2024.png` | Normalized RMSE; shared data-derived color scale |
| Supplement D figures | `replication_{headline,decile,field}_2023.png`, `ais_c_vs_ais_o_six_panel_2023.png` | 2023 common-support replication |
| Supplement D table | `supp_D1_paired.tex` | 2023 paired intervals; annual metrics are also in Main Table 7 |
| Supplement E figures | `maximum_{headline,scatter,decile,field}_2024.png` | Largest treatment-matched samples |
| Supplement E tables | `supp_E1_paired.tex`, `supp_E2_annual.tex` | Maximum-support contrasts and annual metrics |
| Supplement F figures | `leiden_{headline,scatter,decile,field}_2024.png` | Four-cell matched J/L rectangle |
| Supplement F tables | `supp_F1_paired.tex`, `supp_F2_annual.tex` | L robustness contrasts and annual metrics |
| Text: reference gaps by universe/decile | `reference_gaps.csv` | Pooled article counts, not unweighted journal percentages |
| Text: fields below threshold | `field_inventory.csv` | Common-support field counts and n>=300 flag |
| Text: Crossref validation proportions | `crossref_validation.csv` | Classification frequencies in the archived bounded sample |
| Text: coverage loss, isolated/dangling nodes | `coverage.csv` | Counts from supplied full state/node diagnostics |

Functions: `analysis.py` builds supports, descriptors and orchestrates outputs;
`statistics.py` contains the numerical estimators; `contrasts.py` registers the
comparisons and component summaries; `plots.py` and `supplement.py` render figures;
`tables.py` formats the numerical table counterparts.

Source-roster download totals, individual continuation examples, original
matching decisions and metadata-extraction counts are upstream provenance facts.
They require the preserved upstream records; they cannot be derived from final
scores alone. This package does not disguise them as newly reproduced analyses.
