# Replication package for the paper "Reference Observability and Journal Universe Choice in Open Reconstructions of Recursive Journal Metrics" by Utz Weitzel (VU Amsterdam & Radboud University)

**Version 0.1.0 — downstream analysis from precomputed scores.**

This package recreates the statistical analyses, numerical table contents and
figures from journal-level inputs. NF/ANS construction is an upstream step:
the solver is documented and released by [Opindx](https://opindx.org).

## Start with the fictional demonstration

Use Python 3.11 or later. The tested environment is listed in
`requirements-tested.txt`. From this directory:

```console
python -m pip install -r requirements-tested.txt
python -m paper_c_replication demo --work-dir ../paper-c-demo
```

The demo creates its own fictional inputs, then runs the main and supplementary
analyses, including all eight component cells. It uses 199 bootstrap draws for
speed. For a full-sized bootstrap demonstration, add `--draws 4999`.
Figures are marked **SYNTHETIC DEMONSTRATION — NOT PAPER RESULTS**.
The fake journals and values are generated from a fixed seed independently of
empirical observations. They are not perturbed, shuffled, fitted or resampled
Clarivate data. Demo outputs do not reproduce or validate the empirical findings.

## Empirical replication

Prepare the six input CSVs and manifest documented in
[INPUTS.md](docs/INPUTS.md), in a directory **outside this package**. Then run:

```console
python -m paper_c_replication run --inputs ../paper-c-private-inputs --output ../paper-c-private-results
```

Empirical runs use the registered 4,999 paired, field-stratified draws, the
original seed derivation and the updated rolling OpenAlex fields. Inputs are
never discovered automatically, downloaded from Clarivate or transmitted anywhere.
Existing output directories are not overwritten. Expected run time depends on
input size and hardware; bootstrap calculations dominate the empirical run.

**The public package alone cannot reproduce benchmark comparisons.** Empirical
replication requires legitimately held official AIS values, the exact local
matching/roster information and precomputed J/component scores. A Clarivate
subscription alone does not supply the author's resolved crosswalk or frozen
derived scores. See the explicit upstream requirements in `docs/INPUTS.md`.

## What is included

- Data-parameterized statistical and plotting code; no duplicated NF/ANS solver.
- An independent synthetic-data generator for every required input.
- Full precomputed L-universe scores for RAW and REF in both years, using the
  dated CWTS lists and the September OpenAlex inputs. Their population is not
  selected using benchmark availability or JCR membership.
- The additional open N scoring-state nodes omitted by the website's journal
  display filter, so that no N rebuild is required just to restore those nodes.
- Input definitions, a result-to-code map, reproducibility details and release
  checks. N/A scores are obtained from the pinned Opindx release.

## What the archive excludes

No real Clarivate metrics, JCR roster, commercial identifiers, matching tables,
benchmark sample flags, J scores, component scores, ranks, residuals, benchmark
aggregates, empirical plots, manuscript files, notebooks, caches or empirical
validation reports are included. Removing journal names would not make those
derived outputs safe to distribute. This exclusion is intentional even for
files formerly described as “aggregate-only.”

The L data are open-input calculations and contain only the explicitly documented
columns. See [DISCLOSURE_BOUNDARY.md](docs/DISCLOSURE_BOUNDARY.md) for the audit
scope and limitations. Any output from an empirical run stays private: scatter
coordinates, group results and membership-dependent quantities can disclose
restricted information.

## Results and provenance

- [RESULTS_MAP.md](docs/RESULTS_MAP.md): every numerical table and figure.
- [INPUTS.md](docs/INPUTS.md): schemas and replacement instructions.
- [METHODS.md](docs/METHODS.md): samples, rounding, fields and bootstrap settings.
- [VALIDATION.md](docs/VALIDATION.md): verification scope and how to repeat it.
- `code_provenance.json`: source-module provenance and adaptations.
- `data/open/provenance.json`: L data provenance, checks and checksum.

Generated LaTeX files are table fragments requiring `booktabs`. They recreate
numbers and readable headings; final manuscript typesetting remains in the
author's separate TeX project. The package does not edit either manuscript.

## Build the reviewed archive

```console
python -m unittest discover -s tests -v
python tools/audit_open_data.py
python tools/build_release.py
python tools/build_release.py --output ../paper-c-replication-v0.1.0.zip
```

The release builder requires both file names and contents to match
`release_inventory.json`. It fails on added/changed files, links and junctions.
It cannot silently approve a modified package. After any modification, review
the changes and establish a new inventory before releasing. Distribute the
verified ZIP, never the parent research checkout or a working directory.

Code: MIT. Bundled author-created open score tables: CC BY 4.0, attribution to
Utz Weitzel, with the upstream CWTS releases, OpenAlex and the Norwegian Register
(HK-dir; NLOD 2.0) acknowledged.
The fictional example data carry no empirical-data redistribution permission.
