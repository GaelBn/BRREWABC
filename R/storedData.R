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
  reserved <- c("attempt_id", "generation", "job_id", "accepted", "retained", "rejected", "stored_name")
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
    arrow::write_parquet(cbind(metadata, data),
                         file.path(object_dir, paste0(attempt_id, ".parquet")))
  }
  invisible(NULL)
}

schemaSignature <- function(data) {
  paste(paste(names(data), vapply(data, function(x) paste(class(x), collapse = "/"), character(1)), sep = ":"), collapse = "|")
}

persistStoredGeneration <- function(staging_root, storage_root, generation,
                                    retained_ids, summaries_policy, outputs_policy) {
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
      arrow::write_parquet(combined, destination)
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
    arrow::write_parquet(dplyr::bind_rows(previous_manifest, additions), manifest_path)
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

readStoredData <- function(x, kind, name = NULL, generation = NULL,
                           attempt_id = NULL,
                           status = c("all", "accepted", "retained", "rejected")) {
  status <- match.arg(status)
  manifest <- list_abc_stored_data(x)
  if (!nrow(manifest)) return(data.frame())
  selected <- manifest[manifest$kind == kind, , drop = FALSE]
  if (!is.null(name)) selected <- selected[selected$name %in% name, , drop = FALSE]
  if (!is.null(generation)) selected <- selected[selected$generation %in% generation, , drop = FALSE]
  if (!nrow(selected)) return(data.frame())
  tables <- lapply(seq_len(nrow(selected)), function(i) {
    data_path <- selected$file[i]
    if (!grepl("^(/|[A-Za-z]:)", data_path)) data_path <- file.path(resolveStoragePath(x), data_path)
    value <- as.data.frame(arrow::read_parquet(data_path))
    value$stored_name <- selected$name[i]
    value
  })
  result <- dplyr::bind_rows(tables)
  if (!is.null(attempt_id)) result <- result[result$attempt_id %in% attempt_id, , drop = FALSE]
  if (status == "accepted") result <- result[result$accepted, , drop = FALSE]
  if (status == "retained") result <- result[result$retained, , drop = FALSE]
  if (status == "rejected") result <- result[!result$accepted, , drop = FALSE]
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
