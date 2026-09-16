normalizeStoragePolicy <- function(x) {
  match.arg(x, c("none", "retained", "accepted", "all"))
}

normalizeModelResult <- function(x, dist_names) {
  if (is.numeric(x) && is.atomic(x)) {
    distances <- x
    summaries <- list()
    outputs <- list()
  } else if (is.list(x) && "distances" %in% names(x)) {
    distances <- x$distances
    summaries <- x$summaries %||% list()
    outputs <- x$outputs %||% list()
  } else {
    stop("A model must return a numeric distance vector or a named list containing `distances`.", call. = FALSE)
  }

  if (!is.numeric(distances) || length(distances) != length(dist_names) ||
      any(!is.finite(distances))) {
    stop(sprintf("`distances` must contain %d finite numeric value(s).", length(dist_names)), call. = FALSE)
  }
  if (!is.null(names(distances)) && all(nzchar(names(distances)))) {
    if (!setequal(names(distances), dist_names)) {
      stop(sprintf("Named distances must be named: %s.", paste(dist_names, collapse = ", ")), call. = FALSE)
    }
    distances <- distances[dist_names]
  }
  names(distances) <- dist_names

  validateStoredCollection(summaries, "summaries")
  validateStoredCollection(outputs, "outputs")
  list(distances = distances, summaries = summaries, outputs = outputs)
}

`%||%` <- function(x, y) if (is.null(x)) y else x

validateStoredCollection <- function(x, label) {
  if (!is.list(x) || (length(x) && (is.null(names(x)) || any(!nzchar(names(x))) || anyDuplicated(names(x))))) {
    stop(sprintf("`%s` must be a named list with unique names.", label), call. = FALSE)
  }
  if (length(x) && any(!grepl("^[A-Za-z][A-Za-z0-9_.-]*$", names(x)))) {
    stop(sprintf("Names in `%s` may contain only letters, numbers, dots, underscores and hyphens, and must start with a letter.", label), call. = FALSE)
  }
  invisible(TRUE)
}

asStoredTable <- function(x) {
  if (is.matrix(x)) x <- as.data.frame(x, stringsAsFactors = FALSE)
  if (is.data.frame(x)) {
    data <- x
  } else if (is.atomic(x) && is.null(dim(x))) {
    data <- data.frame(
      element_index = seq_along(x),
      element_name = if (is.null(names(x))) rep(NA_character_, length(x)) else names(x),
      value = unname(x),
      stringsAsFactors = FALSE
    )
  } else {
    stop("Stored summaries and outputs must be atomic vectors, matrices, or data frames.", call. = FALSE)
  }
  reserved <- c("attempt_id", "generation", "job_id", "batch_id",
                "attempt_index", "accepted", "retained", "rejected",
                "stored_name")
  if (anyDuplicated(names(data)) || any(names(data) %in% reserved)) {
    stop("A stored table has duplicated or reserved column names.", call. = FALSE)
  }
  if (any(vapply(data, is.list, logical(1)))) {
    stop("List-columns are not supported in stored summaries or outputs.", call. = FALSE)
  }
  data
}

stageStoredCollection <- function(collection, kind, policy, accepted, attempt_id,
                                  generation, job_id, root) {
  if (policy == "none" || (policy %in% c("accepted", "retained") && !accepted) || !length(collection)) {
    return(invisible(NULL))
  }
  for (object_name in names(collection)) {
    data <- asStoredTable(collection[[object_name]])
    if (!nrow(data)) next
    object_dir <- file.path(root, kind, object_name)
    dir.create(object_dir, recursive = TRUE, showWarnings = FALSE)
    metadata <- data.frame(
      attempt_id = rep(attempt_id, nrow(data)),
      generation = rep(as.integer(generation), nrow(data)),
      job_id = rep(as.integer(job_id), nrow(data)),
      accepted = rep(isTRUE(accepted), nrow(data)),
      stringsAsFactors = FALSE
    )
    writeParquetAtomically(
      cbind(metadata, data),
      file.path(object_dir, paste0(attempt_id, ".parquet"))
    )
  }
  invisible(NULL)
}

newStoredBuffer <- function(spec, root, chunk_rows = 1000000L,
                            chunk_mb = 128) {
  chunk_rows <- validatePositiveInteger(chunk_rows, "storage_chunk_rows")
  if (length(chunk_mb) != 1L || is.na(chunk_mb) || !is.finite(chunk_mb) ||
      chunk_mb <= 0) {
    stop("`storage_chunk_mb` must be a positive finite number.",
         call. = FALSE)
  }
  buffer <- new.env(parent = emptyenv())
  buffer$spec <- spec
  buffer$root <- root
  buffer$chunk_rows <- chunk_rows
  buffer$chunk_bytes <- as.double(chunk_mb) * 1024^2
  buffer$entries <- new.env(parent = emptyenv())
  buffer$fragments <- list()
  buffer
}

storedBufferKey <- function(kind, name) paste(kind, name, sep = "\034")

flushStoredBufferEntry <- function(buffer, key) {
  entry <- buffer$entries[[key]]
  if (is.null(entry) || !length(entry$tables)) return(invisible(NULL))
  data <- dplyr::bind_rows(entry$tables)
  part <- entry$next_part
  relative <- file.path(
    entry$kind, entry$name,
    sprintf("generation_%04d", buffer$spec$generation),
    sprintf("batch_%08d", buffer$spec$batch_id),
    sprintf("part_%04d.parquet", part)
  )
  path <- file.path(buffer$root, relative)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  writeParquetAtomically(data, path)
  buffer$fragments[[length(buffer$fragments) + 1L]] <- data.frame(
    kind = entry$kind,
    name = entry$name,
    generation = buffer$spec$generation,
    batch_id = buffer$spec$batch_id,
    part = part,
    file = normalizePath(path, winslash = "/", mustWork = FALSE),
    rows = nrow(data),
    bytes = unname(file.info(path)$size),
    schema = entry$schema,
    stringsAsFactors = FALSE
  )
  entry$tables <- list()
  entry$rows <- 0L
  entry$bytes <- 0
  entry$next_part <- part + 1L
  buffer$entries[[key]] <- entry
  invisible(path)
}

appendStoredCollection <- function(buffer, collection, kind, policy, accepted,
                                   attempt_id, generation, job_id,
                                   attempt_index) {
  if (policy == "none" ||
      (policy %in% c("accepted", "retained") && !accepted) ||
      !length(collection)) {
    return(invisible(buffer))
  }
  for (object_name in names(collection)) {
    data <- asStoredTable(collection[[object_name]])
    if (!nrow(data)) next
    signature <- schemaSignature(data)
    metadata <- data.frame(
      attempt_id = rep(attempt_id, nrow(data)),
      generation = rep(as.integer(generation), nrow(data)),
      job_id = rep(as.integer(job_id), nrow(data)),
      batch_id = rep(as.integer(buffer$spec$batch_id), nrow(data)),
      attempt_index = rep(as.integer(attempt_index), nrow(data)),
      accepted = rep(isTRUE(accepted), nrow(data)),
      stringsAsFactors = FALSE
    )
    table <- cbind(metadata, data)
    key <- storedBufferKey(kind, object_name)
    entry <- buffer$entries[[key]]
    if (is.null(entry)) {
      entry <- list(
        kind = kind, name = object_name, schema = signature,
        tables = list(), rows = 0L, bytes = 0, next_part = 1L
      )
    } else if (!identical(entry$schema, signature)) {
      stop(sprintf(
        "The schema of `%s` `%s` changed within batch %d.",
        kind, object_name, buffer$spec$batch_id
      ), call. = FALSE)
    }
    entry$tables[[length(entry$tables) + 1L]] <- table
    entry$rows <- entry$rows + nrow(table)
    entry$bytes <- entry$bytes + as.numeric(utils::object.size(table))
    buffer$entries[[key]] <- entry
    if (entry$rows >= buffer$chunk_rows ||
        entry$bytes >= buffer$chunk_bytes) {
      flushStoredBufferEntry(buffer, key)
    }
  }
  invisible(buffer)
}

flushStoredBuffer <- function(buffer) {
  keys <- ls(buffer$entries, all.names = TRUE)
  for (key in keys) flushStoredBufferEntry(buffer, key)
  dplyr::bind_rows(buffer$fragments)
}

schemaSignature <- function(data) {
  paste(paste(names(data), vapply(data, function(x) paste(class(x), collapse = "/"), character(1)), sep = ":"), collapse = "|")
}

writeParquetAtomically <- function(data, path) {
  temporary_path <- paste0(path, ".tmp-", Sys.getpid())
  on.exit(unlink(temporary_path), add = TRUE)
  arrow::write_parquet(data, temporary_path)
  if (!file.rename(temporary_path, path)) {
    stop(sprintf("Unable to finalize Parquet file `%s`.", path), call. = FALSE)
  }
  invisible(path)
}

manifestActiveRows <- function(manifest) {
  if (is.null(manifest) || !nrow(manifest)) return(manifest)
  if ("active" %in% names(manifest)) {
    manifest[is.na(manifest$active) | manifest$active, , drop = FALSE]
  } else {
    manifest
  }
}

persistBatchFragments <- function(fragments, storage_root, generation,
                                  committed_ids, retained_ids,
                                  summaries_policy, outputs_policy) {
  manifest_path <- file.path(storage_root, "manifest.parquet")
  previous_manifest <- if (file.exists(manifest_path)) {
    as.data.frame(arrow::read_parquet(manifest_path))
  } else {
    NULL
  }
  active_manifest <- manifestActiveRows(previous_manifest)
  additions <- list()

  if (is.null(fragments) || !nrow(fragments)) return(invisible(data.frame()))
  fragments <- fragments[fragments$generation == generation, , drop = FALSE]
  for (index in seq_len(nrow(fragments))) {
    fragment <- fragments[index, , drop = FALSE]
    policy <- if (fragment$kind == "summaries") {
      summaries_policy
    } else {
      outputs_policy
    }
    if (policy == "none" || !file.exists(fragment$file)) next
    data <- as.data.frame(arrow::read_parquet(fragment$file))
    data <- data[data$attempt_id %in% committed_ids, , drop = FALSE]
    if (policy == "retained") {
      data <- data[data$attempt_id %in% retained_ids, , drop = FALSE]
    }
    if (!nrow(data)) next

    metadata_names <- c(
      "attempt_id", "generation", "job_id", "batch_id",
      "attempt_index", "accepted", "retained"
    )
    signature <- schemaSignature(data[setdiff(names(data), metadata_names)])
    old_schema <- if (!is.null(active_manifest)) {
      unique(active_manifest$schema[
        active_manifest$kind == fragment$kind &
          active_manifest$name == fragment$name
      ])
    } else {
      character()
    }
    new_schema <- vapply(
      additions,
      function(x) if (x$kind == fragment$kind && x$name == fragment$name) {
        x$schema
      } else {
        NA_character_
      },
      character(1)
    )
    expected_schema <- unique(c(old_schema, stats::na.omit(new_schema)))
    if (!identical(signature, fragment$schema) ||
        (length(expected_schema) && !identical(signature, expected_schema))) {
      stop(sprintf(
        "The schema of `%s` `%s` is not fixed across simulations.",
        fragment$kind, fragment$name
      ), call. = FALSE)
    }

    data$retained <- data$attempt_id %in% retained_ids
    ordered_metadata <- intersect(metadata_names, names(data))
    data <- data[c(
      ordered_metadata,
      setdiff(names(data), ordered_metadata)
    )]
    destination_relative <- file.path(
      fragment$kind, fragment$name,
      sprintf("generation_%04d", generation),
      sprintf(
        "part-batch%08d-%04d.parquet",
        fragment$batch_id, fragment$part
      )
    )
    destination <- file.path(storage_root, destination_relative)
    dir.create(dirname(destination), recursive = TRUE, showWarnings = FALSE)
    writeParquetAtomically(data, destination)
    additions[[length(additions) + 1L]] <- data.frame(
      format_version = 2L,
      kind = fragment$kind,
      name = fragment$name,
      generation = as.integer(generation),
      layout = "fragmented",
      batch_id = as.integer(fragment$batch_id),
      part = as.integer(fragment$part),
      file = destination_relative,
      rows = nrow(data),
      bytes = unname(file.info(destination)$size),
      schema = signature,
      active = TRUE,
      stringsAsFactors = FALSE
    )
  }

  additions <- dplyr::bind_rows(additions)
  if (!nrow(additions)) return(invisible(additions))
  if (!is.null(previous_manifest)) {
    replaced_keys <- unique(paste(
      additions$kind, additions$name, additions$generation
    ))
    previous_manifest <- previous_manifest[
      !paste(
        previous_manifest$kind, previous_manifest$name,
        previous_manifest$generation
      ) %in% replaced_keys,
      , drop = FALSE
    ]
  }
  writeParquetAtomically(
    dplyr::bind_rows(previous_manifest, additions), manifest_path
  )
  unlink(unique(fragments$file))
  invisible(additions)
}

persistStoredGeneration <- function(staging_root, storage_root, generation,
                                    committed_ids, retained_ids,
                                    summaries_policy, outputs_policy,
                                    fragments = NULL) {
  if (!is.null(fragments)) {
    return(persistBatchFragments(
      fragments, storage_root, generation, committed_ids, retained_ids,
      summaries_policy, outputs_policy
    ))
  }
  manifest_path <- file.path(storage_root, "manifest.parquet")
  previous_manifest <- if (file.exists(manifest_path)) as.data.frame(arrow::read_parquet(manifest_path)) else NULL
  new_manifest <- list()

  for (kind in c("summaries", "outputs")) {
    policy <- if (kind == "summaries") summaries_policy else outputs_policy
    kind_root <- file.path(staging_root, kind)
    if (policy == "none" || !dir.exists(kind_root)) next
    object_dirs <- list.dirs(kind_root, recursive = FALSE, full.names = TRUE)
    for (object_dir in object_dirs) {
      files <- list.files(object_dir, pattern = "\\.parquet$", full.names = TRUE)
      if (!length(files)) next
      records <- lapply(files, function(path) as.data.frame(arrow::read_parquet(path)))
      records <- records[vapply(records, function(z) all(z$generation == generation), logical(1))]
      records <- records[vapply(records, function(z) z$attempt_id[1] %in% committed_ids, logical(1))]
      if (!length(records)) next
      if (policy == "retained") records <- records[vapply(records, function(z) z$attempt_id[1] %in% retained_ids, logical(1))]
      if (!length(records)) next

      metadata_names <- c("attempt_id", "generation", "job_id", "accepted")
      signatures <- unique(vapply(records, function(z) schemaSignature(z[setdiff(names(z), metadata_names)]), character(1)))
      object_name <- basename(object_dir)
      old_schema <- if (!is.null(previous_manifest)) unique(previous_manifest$schema[previous_manifest$kind == kind & previous_manifest$name == object_name]) else character()
      if (length(signatures) != 1L || (length(old_schema) && !identical(signatures, old_schema))) {
        stop(sprintf("The schema of `%s` `%s` is not fixed across simulations.", kind, object_name), call. = FALSE)
      }

      tables <- lapply(records, function(z) {
        z$retained <- z$attempt_id %in% retained_ids
        z[c("attempt_id", "generation", "job_id", "accepted", "retained",
            setdiff(names(z), c(metadata_names, "retained")))]
      })
      combined <- dplyr::bind_rows(tables)
      destination_relative <- file.path(kind, object_name, sprintf("generation_%04d.parquet", generation))
      destination_dir <- file.path(storage_root, kind, object_name)
      dir.create(destination_dir, recursive = TRUE, showWarnings = FALSE)
      destination <- file.path(storage_root, destination_relative)
      writeParquetAtomically(combined, destination)
      new_manifest[[length(new_manifest) + 1L]] <- data.frame(
        kind = kind, name = object_name, generation = as.integer(generation),
        file = destination_relative,
        rows = nrow(combined), schema = signatures, stringsAsFactors = FALSE
      )
    }
  }

  additions <- dplyr::bind_rows(new_manifest)
  if (nrow(additions)) {
    dir.create(storage_root, recursive = TRUE, showWarnings = FALSE)
    if (!is.null(previous_manifest)) {
      keys <- paste(additions$kind, additions$name, additions$generation)
      previous_manifest <- previous_manifest[!paste(previous_manifest$kind, previous_manifest$name, previous_manifest$generation) %in% keys, , drop = FALSE]
    }
    writeParquetAtomically(dplyr::bind_rows(previous_manifest, additions), manifest_path)
  }
  invisible(additions)
}

resolveStoragePath <- function(x) {
  if (is.list(x) && !is.null(x$storage$path)) return(x$storage$path)
  if (!is.character(x) || length(x) != 1L) stop("`x` must be an ABC result or an experiment directory.", call. = FALSE)
  direct <- file.path(x, "manifest.parquet")
  nested <- file.path(x, "res", "parquet", "manifest.parquet")
  if (file.exists(direct)) x else dirname(nested)
}

#' List summary statistics and outputs stored for an ABC result
#'
#' @param x an object returned by [abcsmc()] or [abcrejection()], or the path to
#' an experiment directory or its `res/parquet` directory.
#' @return a data frame describing available Parquet datasets.
#' @export
list_abc_stored_data <- function(x) {
  root <- resolveStoragePath(x)
  manifest <- file.path(root, "manifest.parquet")
  if (!file.exists(manifest)) return(data.frame())
  as.data.frame(arrow::read_parquet(manifest))
}

storageFilePaths <- function(root, files) {
  vapply(files, function(path) {
    if (grepl("^(/|[A-Za-z]:)", path)) path else file.path(root, path)
  }, character(1))
}

storageRelativePath <- function(root, path) {
  root <- normalizePath(root, winslash = "/", mustWork = FALSE)
  path <- normalizePath(path, winslash = "/", mustWork = FALSE)
  prefix <- paste0(root, "/")
  if (startsWith(path, prefix)) substring(path, nchar(prefix) + 1L) else path
}

writeParquetStreamAtomically <- function(paths, destination,
                                         compression = "zstd",
                                         chunk_size = 1048576L,
                                         overwrite = FALSE) {
  if (!length(paths)) stop("No Parquet fragments to consolidate.", call. = FALSE)
  dir.create(dirname(destination), recursive = TRUE, showWarnings = FALSE)
  if (file.exists(destination) && !overwrite) {
    stop(sprintf("File `%s` already exists.", destination), call. = FALSE)
  }
  temporary <- tempfile(
    pattern = paste0(basename(destination), ".tmp-"),
    tmpdir = dirname(destination)
  )
  on.exit(unlink(temporary), add = TRUE)
  first <- arrow::read_parquet(paths[[1L]], as_data_frame = FALSE)
  sink <- arrow::FileOutputStream$create(temporary)
  writer <- arrow::ParquetFileWriter$create(
    first$schema,
    sink,
    arrow::ParquetWriterProperties$create(
      column_names = first$ColumnNames(), compression = compression
    )
  )
  writer_open <- TRUE
  sink_open <- TRUE
  on.exit({
    if (writer_open) try(writer$Close(), silent = TRUE)
    if (sink_open) try(sink$close(), silent = TRUE)
  }, add = TRUE)
  for (path in paths) {
    table <- arrow::read_parquet(path, as_data_frame = FALSE)
    writer$WriteTable(table, as.integer(chunk_size))
  }
  writer$Close()
  writer_open <- FALSE
  sink$close()
  sink_open <- FALSE
  if (file.exists(destination)) unlink(destination)
  if (!file.rename(temporary, destination)) {
    stop(sprintf("Unable to finalize file `%s`.", destination),
         call. = FALSE)
  }
  invisible(destination)
}

#' Consolidate stored ABC summaries and outputs
#'
#' Consolidates the active Parquet fragments for each selected named object and
#' generation into one Parquet file. By default this creates a separate export
#' and leaves the operational fragmented storage unchanged.
#'
#' @param x an object returned by [abcsmc()] or [abcrejection()], or an
#' experiment/storage directory.
#' @param kind one or both of `"summaries"` and `"outputs"`.
#' @param name optional stored object names.
#' @param generation optional generation numbers.
#' @param mode `"export"` to leave the active manifest unchanged, or
#' `"replace"` to make consolidated files the active storage layout.
#' @param output_dir destination directory. Defaults to `consolidated` for an
#' export and `compacted` for replacement, below the storage root.
#' @param keep_fragments whether source fragments should remain on disk after a
#' successful replacement.
#' @param compression Parquet compression codec.
#' @param overwrite whether existing consolidated files may be replaced.
#' @param dry_run return the consolidation plan without writing files.
#' @return a data frame describing the consolidation performed or planned.
#' @export
consolidate_abc_storage <- function(
    x,
    kind = c("summaries", "outputs"),
    name = NULL,
    generation = NULL,
    mode = c("export", "replace"),
    output_dir = NULL,
    keep_fragments = TRUE,
    compression = "zstd",
    overwrite = FALSE,
    dry_run = FALSE) {
  mode <- match.arg(mode)
  kind <- match.arg(kind, c("summaries", "outputs"), several.ok = TRUE)
  root <- resolveStoragePath(x)
  manifest_path <- file.path(root, "manifest.parquet")
  manifest <- list_abc_stored_data(x)
  if (!nrow(manifest)) return(data.frame())
  selected <- manifestActiveRows(manifest)
  selected <- selected[selected$kind %in% kind, , drop = FALSE]
  if (!is.null(name)) selected <- selected[selected$name %in% name, , drop = FALSE]
  if (!is.null(generation)) {
    selected <- selected[selected$generation %in% generation, , drop = FALSE]
  }
  if (!nrow(selected)) return(data.frame())
  if (is.null(output_dir)) {
    output_dir <- file.path(
      root, if (mode == "export") "consolidated" else "compacted"
    )
  }
  output_dir <- normalizePath(output_dir, winslash = "/", mustWork = FALSE)
  group_key <- paste(selected$kind, selected$name, selected$generation, sep = "\034")
  groups <- split(seq_len(nrow(selected)), group_key)
  plan <- lapply(groups, function(indices) {
    rows <- selected[indices, , drop = FALSE]
    paths <- storageFilePaths(root, rows$file)
    destination <- file.path(
      output_dir, rows$kind[[1L]], rows$name[[1L]],
      sprintf("generation_%04d.parquet", rows$generation[[1L]])
    )
    data.frame(
      kind = rows$kind[[1L]],
      name = rows$name[[1L]],
      generation = as.integer(rows$generation[[1L]]),
      source_files = length(paths),
      source_rows = sum(rows$rows),
      source_bytes = sum(file.info(paths)$size, na.rm = TRUE),
      output_file = normalizePath(
        destination, winslash = "/", mustWork = FALSE
      ),
      stringsAsFactors = FALSE
    )
  })
  plan <- dplyr::bind_rows(plan)
  if (dry_run) return(plan)

  plan$output_bytes <- numeric(nrow(plan))
  plan$replaced <- logical(nrow(plan))
  additions <- list()
  source_paths <- character()
  for (index in seq_len(nrow(plan))) {
    matching <- selected$kind == plan$kind[[index]] &
      selected$name == plan$name[[index]] &
      selected$generation == plan$generation[[index]]
    rows <- selected[matching, , drop = FALSE]
    paths <- storageFilePaths(root, rows$file)
    signatures <- unique(rows$schema)
    if (length(signatures) != 1L) {
      stop(sprintf(
        "Cannot consolidate `%s` `%s`: incompatible schemas.",
        plan$kind[[index]], plan$name[[index]]
      ), call. = FALSE)
    }
    writeParquetStreamAtomically(
      paths, plan$output_file[[index]], compression = compression,
      overwrite = overwrite
    )
    plan$output_bytes[index] <- unname(
      file.info(plan$output_file[[index]])$size
    )
    plan$replaced[index] <- identical(mode, "replace")
    if (mode == "replace") {
      additions[[length(additions) + 1L]] <- data.frame(
        format_version = 2L,
        kind = plan$kind[[index]],
        name = plan$name[[index]],
        generation = plan$generation[[index]],
        layout = "consolidated",
        batch_id = NA_integer_,
        part = NA_integer_,
        file = storageRelativePath(root, plan$output_file[[index]]),
        rows = plan$source_rows[[index]],
        bytes = plan$output_bytes[[index]],
        schema = signatures,
        active = TRUE,
        stringsAsFactors = FALSE
      )
      source_paths <- c(source_paths, paths)
    }
  }
  if (mode == "replace") {
    selected_keys <- unique(paste(
      selected$kind, selected$name, selected$generation
    ))
    remaining <- manifest[
      !paste(manifest$kind, manifest$name, manifest$generation) %in%
        selected_keys,
      , drop = FALSE
    ]
    writeParquetAtomically(
      dplyr::bind_rows(remaining, dplyr::bind_rows(additions)),
      manifest_path
    )
    if (!keep_fragments) {
      destinations <- normalizePath(
        plan$output_file, winslash = "/", mustWork = FALSE
      )
      sources <- normalizePath(
        unique(source_paths), winslash = "/", mustWork = FALSE
      )
      unlink(setdiff(sources, destinations))
    }
  }
  plan
}

readStoredData <- function(x, kind, name = NULL, generation = NULL,
                           attempt_id = NULL,
                           status = c("all", "accepted", "retained", "rejected")) {
  status <- match.arg(status)
  manifest <- list_abc_stored_data(x)
  if (!nrow(manifest)) return(data.frame())
  selected <- manifestActiveRows(manifest)
  selected <- selected[selected$kind == kind, , drop = FALSE]
  if (!is.null(name)) selected <- selected[selected$name %in% name, , drop = FALSE]
  if (!is.null(generation)) selected <- selected[selected$generation %in% generation, , drop = FALSE]
  if (!nrow(selected)) return(data.frame())
  root <- resolveStoragePath(x)
  tables <- lapply(unique(selected$name), function(object_name) {
    object_rows <- selected[selected$name == object_name, , drop = FALSE]
    paths <- storageFilePaths(root, object_rows$file)
    query <- arrow::open_dataset(paths, format = "parquet")
    if (!is.null(attempt_id)) {
      requested_ids <- attempt_id
      query <- dplyr::filter(query, attempt_id %in% !!requested_ids)
    }
    if (status == "accepted") query <- dplyr::filter(query, accepted)
    if (status == "retained") query <- dplyr::filter(query, retained)
    if (status == "rejected") query <- dplyr::filter(query, !accepted)
    value <- as.data.frame(dplyr::collect(query))
    value$stored_name <- rep(object_name, nrow(value))
    value
  })
  result <- dplyr::bind_rows(tables)
  rownames(result) <- NULL
  result
}

#' Read summary statistics stored in Parquet files
#'
#' @inheritParams list_abc_stored_data
#' @param name optional names of summary statistics to read.
#' @param generation optional generation numbers.
#' @param attempt_id optional attempt identifiers.
#' @param status one of `"all"`, `"accepted"`, `"retained"`, or `"rejected"`.
#' @return a data frame containing the selected summary statistics.
#' @export
read_summary_statistics <- function(x, name = NULL, generation = NULL,
                                    attempt_id = NULL,
                                    status = c("all", "accepted", "retained", "rejected")) {
  readStoredData(x, "summaries", name, generation, attempt_id, status)
}

#' Read model outputs stored in Parquet files
#'
#' @inheritParams read_summary_statistics
#' @return a data frame containing the selected model outputs.
#' @export
read_model_outputs <- function(x, name = NULL, generation = NULL,
                               attempt_id = NULL,
                               status = c("all", "accepted", "retained", "rejected")) {
  readStoredData(x, "outputs", name, generation, attempt_id, status)
}

utils::globalVariables(c("accepted", "retained"))
