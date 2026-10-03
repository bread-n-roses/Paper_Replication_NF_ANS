# Distribution boundary

The package follows the author's explicit instruction to exclude Clarivate
data and information that could reconstruct them. It therefore treats real
J-derived outputs as private even when the values were calculated from open
citations. A list of J scores would itself disclose the roster. A residual plus
an open score can reveal the benchmark. Renaming identifiers, hashing journal
IDs, rounding, or calling a table “aggregate-only” is not our release strategy.

## Permitted contents

- Reviewed source and documentation describing algorithms and input schemas.
- A deterministic synthetic generator that does not read or fit empirical data.
- Full L score tables derived from OpenAlex and publicly released CWTS lists,
  without any filter or join based on JCR/Clarivate.
- Additional N nodes omitted by Opindx's journal display filter, selected solely
  by that public-export difference, without JCR/Clarivate conditions.
- Checksums of public package files and source-code provenance.

## Excluded contents

All actual official metrics; JCR lists and identifiers; private crosswalks;
J state membership/scores; component scores; comparison sample membership;
official-score missingness flags, ranks and deciles; errors and residuals;
benchmark-conditioned subsets of open data; actual reference-validation sample
rows; benchmark aggregates; plots, TeX/PDF manuscripts, logs, caches, pickles,
notebooks, empirical validation files and hashes of restricted input files.

Empirical graphics are excluded even if rasterized and anonymous. Coordinates
may be recoverable, and combining results can enable further inference. Public
L data are the entire L scoring population; their inclusion does not identify
which journals were in any benchmark comparison.

## Controls

1. Default-deny package `.gitignore`; explicit source/data file selection.
2. Release inventory pins exact bytes, not just filename extensions.
3. Builder rejects unknown files, changed hashes, path escapes, links and junctions.
4. Work directories must resolve outside the package directory. Empirical outputs
   are always marked `DO_NOT_DISTRIBUTE`.
5. No local configuration, parent-directory discovery, executable notebook output,
   credentials, or empirical values are required for the fictional demo.
6. The final ZIP is opened again and every entry is compared to the inventory.
7. Private author-side checks compare runtime results with the accepted analysis;
   the detailed checks and expected benchmark values are kept outside the archive.

These controls apply to the reviewed archive, not arbitrary future working
directories. The inventory is a reproducibility/integrity control, not a digital
signature or an automatic confidentiality classifier. Editing both source and
inventory requires a new content review. A future user can still put confidential
information in a permitted-looking source file; a filename check alone cannot
establish safety.

## Limits of the assurance

The technical claim is that the reviewed distribution contains no empirical
Clarivate input or benchmark-conditioned result. It provides no transformation
or encoded copy from which those inputs can be inverted. Open scores may be
statistically associated with other metrics; an assertion of zero possible
inference under all external information would be stronger than this audit can
establish. This is not a legal determination of rights under an institutional
contract, and it does not certify the separately maintained paper or supplement.

QSS explicitly allows justified proprietary-data exceptions in its
[submission guidance](https://direct.mit.edu/qss/pages/submission-guidelines).
The package makes its restricted-input dependency explicit; it does not infer
redistribution permission from that journal policy. Clarivate's own
[JCR usage guidance](https://clarivate.com/academia-government/blog/a-quick-refresher-on-journal-citation-reports-use-cases-branding-and-terms-of-use/)
is separate from any applicable institutional agreement. The package avoids
relying on permission to redistribute empirical Clarivate-derived results.
