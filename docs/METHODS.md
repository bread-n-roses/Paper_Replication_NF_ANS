# Registered analytical settings

- Years: 2023 and 2024; primary comparisons use 2024.
- RAW/REF and J/N/A form six cells. L is a separate robustness universe.
- Common support: nonmissing official AIS and nonmissing ANS in all six cells.
- Maximum treatment support: official AIS and both treatments in the specified
  universe. Maximum universe support: official AIS and the two named cells.
- L plots use the four-cell J/L rectangle; paired contrast intervals use their
  stated maximum pairwise supports. These are not interchangeable.
- Undefined values are never replaced by zero. Genuine numerical zeroes remain.
- MAE = mean absolute ANS–AIS difference; RMSE = square root of mean squared
  difference. Signed error = mean ANS minus AIS. Spearman uses average ranks.
- Scaled MAE and NRMSE divide group error by mean official AIS in that same group.
- Official AIS deciles are recomputed on the stated comparison sample using
  maximum-rank percentiles and preserving ties. D1 is highest. Top decile is
  strictly above the 90th percentile.
- Field is the score-year rolling modal OpenAlex primary-work field. Missing
  fields become `Unknown`. Display fields with at least 300 journals, computed
  on the applicable sample; excluded fields remain in machine-readable outputs.
- Paired changes compare scores for the same journals. Mean absolute score
  movement is different from the change in benchmark MAE.
- Bootstrap: 4,999 paired journal resamples, within field, retaining each field's
  observed sample size. Whole rows remain paired. The base seed is 20260821;
  each scalar gets its own SHA256-derived child seed from
  `year|support_id|contrast_id|statistic`. Quantiles 0.025 and 0.975 give intervals.
- The package computes the 48 scalar intervals printed in the accepted main
  paper and supplement. Other contrast point estimates are also produced. It
  does not claim to rerun unreported legacy 2025, EF-bridge or rank-bootstrap
  diagnostics.
- Component effects: journal-level factorial contrasts across eight supplied
  J-universe ANS values; average absolute contrasts only after calculating them
  journal by journal. Internal A = manuscript D (denominator).
- Full precision is retained throughout. Final numerical tables normally display
  three decimals; heatmap annotations display two. Counts and displayed
  percentages may need the manuscript's corresponding formatting.

The statistical module preserves the original computation and resampling order.
The package does not tune synthetic distributions to empirical results, choose
universes based on benchmark agreement, or recompute upstream journal metrics.
