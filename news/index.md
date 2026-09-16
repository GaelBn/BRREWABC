# Changelog

## BRREWABC 1.3.0

- Added optional storage of model summary statistics and detailed
  outputs for ABC rejection and ABC-SMC analyses.
- Model functions can now return a structured list containing
  `distances`, `summaries`, and `outputs`; the historical numeric
  distance vector remains supported.
- Added Parquet storage for vectors, matrices, and data frames, with
  independent `none`, `retained`, `accepted`, and `all` retention
  policies for summaries and outputs.
- Added stable attempt identifiers and acceptance/retention metadata to
  link stored tables to tested particles.
- Parquet fragments are finalized atomically so stopping parallel
  workers cannot expose partially written files during generation
  consolidation.
- A tested particle is now committed before it is published as accepted,
  and Parquet consolidation ignores orphan fragments not present in the
  committed attempt journal.
- Added
  [`list_abc_stored_data()`](https://gaelbn.github.io/BRREWABC/reference/list_abc_stored_data.md),
  [`read_summary_statistics()`](https://gaelbn.github.io/BRREWABC/reference/read_summary_statistics.md),
  and
  [`read_model_outputs()`](https://gaelbn.github.io/BRREWABC/reference/read_model_outputs.md)
  to inspect and selectively load stored data by name, generation,
  attempt identifier, or ABC status.

## BRREWABC 1.2.1

- Fixed a bug when many plots are generated.

## BRREWABC 1.2.0

- Update parallel task schedulling

## BRREWABC 1.1.0

- Add ABC rejection method and a function to compute ESS (ABC-SMC)

## BRREWABC 1.0.0

- First public release.
