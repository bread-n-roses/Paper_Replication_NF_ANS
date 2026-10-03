# Validation and its limits

## Verification performed for version 0.1.0

The author ran the portable analysis locally on the frozen empirical inputs and
compared the results with the accepted September/rolling-field analysis used in
the version-4 manuscript and supplement. The private comparison passed 23
groups of checks, including:

- Primary annual statistics, coverage and sample composition.
- Reference-observability counts and all named contrast point estimates.
- All 48 printed confidence intervals, with the registered 4,999 draws.
- Component-cell benchmark statistics and journal-level factorial effects.
- Decile profiles and the accepted rolling-field profiles.
- Maximum-sample and Leiden robustness statistics and field profiles.
- Exact comparison-sample membership for both years, not only sample sizes.

The detailed report, expected empirical values, hashes of private inputs, and
all empirical results remain outside the public package. This is an author-side
verification statement; a reader without the restricted inputs cannot
independently validate those empirical comparisons from the public demo.

The fictional workflow was run separately and generated the 18 main/supplementary
figures, machine-readable numerical summaries and LaTeX table fragments. The
demo has no access to real benchmark values. Representative scatter and field
plots were visually inspected; this is not a claim of pixel-identical manuscript
typesetting. LaTeX fragments are not compiled into the author's paper here.

## Repeat the public checks

```console
python -m unittest discover -s tests -v
python tools/audit_open_data.py
python tools/build_release.py
python -m paper_c_replication demo --work-dir ../fresh-demo
```

Tests cover known factorial contrasts, the distinction between benchmark error
and score movement, rank/decile ties, paired-bootstrap reproducibility,
exact CSV round trips, duplicate input rejection, no empirical data reads by
the synthetic generator, output-directory isolation and release rejection of
changed/extra files. The open-data audit validates exact columns, checksums,
complete RAW/REF state agreement, count ordering and NF/ANS identities.

## Empirical checks for licensed users

Run the empirical command with the exact six input files described in INPUTS.md.
Compare numerical outputs at full precision to your preserved accepted-run
outputs, and compare displayed values at the manuscript's rounding precision.
Check support membership and definitions before interpreting any numerical
discrepancy. Match the pinned input versions, field labels, random seed and
bootstrap draw count. CSV loading uses `float_precision='round_trip'` to avoid
changing exact score ties.

A successful run creates `run_report.json` with software versions and local
input hashes. The flag `empirical_results_validated` remains false: execution
alone is not a comparison against an external empirical reference. Only a
separate comparison justifies that conclusion. Keep the report and all empirical
outputs private; the release builder does not include them.
