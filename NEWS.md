# BRREWABC 1.3.0

* Replaced concurrent appends to shared CSV files and their file locks with a
  coordinator-owned, batched execution protocol. Workers now write private
  batch results atomically and local batches run in isolated `callr` processes.
* Added `batch_size` to control how many simulations each worker performs before
  returning to the coordinator.
* Summary statistics and model outputs are now buffered and written as bounded,
  atomic Parquet fragments per batch instead of one temporary file per particle.
  Fragment metadata is committed with the batch result, and generation
  finalization uses bounded memory without rescanning previous generations.
* Added `storage_chunk_rows` and `storage_chunk_mb` to bound worker storage
  buffers independently of the simulation batch size.
* Added `consolidate_abc_storage()` to export or transactionally replace active
  fragments with one Parquet file per named object and generation.
* Added optional storage of model summary statistics and detailed outputs for
  ABC rejection and ABC-SMC analyses.
* Model functions can now return a structured list containing `distances`,
  `summaries`, and `outputs`; the historical numeric distance vector remains
  supported.
* Added Parquet storage for vectors, matrices, and data frames, with independent
  `none`, `retained`, `accepted`, and `all` retention policies for summaries and
  outputs.
* Added stable attempt identifiers and acceptance/retention metadata to link
  stored tables to tested particles.
* Parquet fragments are finalized atomically so stopping parallel workers cannot
  expose partially written files.
* A tested particle is now committed before it is published as accepted, and
  Parquet consolidation ignores orphan fragments not present in the committed
  attempt journal.
* Added `list_abc_stored_data()`, `read_summary_statistics()`, and
  `read_model_outputs()` to inspect and selectively load stored data by name,
  generation, attempt identifier, or ABC status.

# BRREWABC 1.2.1

* Fixed a bug when many plots are generated.

# BRREWABC 1.2.0

* Update parallel task schedulling

# BRREWABC 1.1.0

* Add ABC rejection method and a function to compute ESS (ABC-SMC)

# BRREWABC 1.0.0

* First public release.
