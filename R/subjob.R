writeRDSAtomically <- function(value, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary_path <- tempfile(
    pattern = paste0(basename(path), ".tmp-"),
    tmpdir = dirname(path)
  )
  on.exit(unlink(temporary_path), add = TRUE)
  saveRDS(value, temporary_path)
  if (!file.rename(temporary_path, path)) {
    stop(sprintf("Unable to finalize file `%s`.", path), call. = FALSE)
  }
  invisible(path)
}

validateBatchSize <- function(batch_size) {
  if (length(batch_size) != 1L || is.na(batch_size) ||
      batch_size < 1 || batch_size != as.integer(batch_size)) {
    stop("`batch_size` must be a single positive integer.", call. = FALSE)
  }
  as.integer(batch_size)
}

validatePositiveInteger <- function(value, name) {
  if (length(value) != 1L || is.na(value) || value < 1 ||
      value != as.integer(value)) {
    stop(sprintf("`%s` must be a positive integer.", name), call. = FALSE)
  }
  as.integer(value)
}

newBatchSpec <- function(batch_id, generation, job_id, n_attempts,
                         lhs_indices = NULL, seed = NULL) {
  if (is.null(seed)) seed <- sample.int(.Machine$integer.max, 1L)
  list(
    batch_id = as.integer(batch_id),
    generation = as.integer(generation),
    job_id = as.integer(job_id),
    n_attempts = as.integer(n_attempts),
    lhs_indices = lhs_indices,
    seed = as.integer(seed)
  )
}

batchAttemptId <- function(spec, attempt_index) {
  sprintf(
    "g%04d-b%08d-a%04d",
    spec$generation, spec$batch_id, attempt_index
  )
}

buildParticleRow <- function(values, metadata, distances, column_names) {
  row <- data.frame(c(metadata, values, distances), check.names = FALSE)
  missing_columns <- setdiff(column_names, colnames(row))
  row[missing_columns] <- NA
  row[column_names]
}

runSMCBatch <- function(spec, state) {
  rows <- vector("list", spec$n_attempts)

  for (attempt_index in seq_len(spec$n_attempts)) {
    lhs_index <- if (length(spec$lhs_indices)) {
      spec$lhs_indices[[attempt_index]]
    } else {
      attempt_index
    }
    attempt_id <- batchAttemptId(spec, attempt_index)
    proposed_particle <- createParticleForSMC(
      state$gen, lhs_index, state$lhs_first_gen,
      state$previous_acc_particles, state$empirical_sd,
      state$use_lhs_for_first_iter, state$prior_dist,
      state$model_names, state$model_jump_prob
    )
    proposed_particle[["job_id"]] <- spec$job_id
    model_result <- normalizeModelResult(
      state$model_list[[proposed_particle[["model"]]]](
        proposed_particle, state$ss_obs
      ),
      state$dist_names
    )
    distances <- model_result$distances
    accepted <- (state$gen == 1L) || all(distances < state$epsilon)

    stageStoredCollection(
      model_result$summaries, "summaries", state$store_summaries,
      accepted, attempt_id, state$gen, spec$job_id,
      state$tmp_object_store_root
    )
    stageStoredCollection(
      model_result$outputs, "outputs", state$store_outputs,
      accepted, attempt_id, state$gen, spec$job_id,
      state$tmp_object_store_root
    )

    p_weight <- NA_real_
    if (state$gen == 1L) {
      p_weight <- 1
    } else if (accepted) {
      p_weight <- computeWeight(
        proposed_particle, state$previous_acc_particles,
        state$empirical_sd, state$prior_dist
      )
    }

    particle_values <- proposed_particle[
      setdiff(names(proposed_particle), "job_id")
    ]
    rows[[attempt_index]] <- buildParticleRow(
      particle_values,
      list(
        gen = state$gen,
        attempt_id = attempt_id,
        job_id = spec$job_id,
        accepted = accepted,
        retained = FALSE,
        pWeight = p_weight,
        batch_id = spec$batch_id,
        attempt_index = attempt_index
      ),
      stats::setNames(as.list(distances), state$dist_names),
      c(state$column_names, "batch_id", "attempt_index")
    )
  }

  dplyr::bind_rows(rows)
}

runRejectionBatch <- function(spec, state) {
  rows <- vector("list", spec$n_attempts)

  for (attempt_index in seq_len(spec$n_attempts)) {
    attempt_id <- batchAttemptId(spec, attempt_index)
    proposed_particle <- createParticle(state$prior_dist, state$model_names)
    proposed_particle[["job_id"]] <- spec$job_id
    model_result <- normalizeModelResult(
      state$model_list[[proposed_particle[["model"]]]](
        proposed_particle, state$ss_obs
      ),
      state$dist_names
    )
    distances <- model_result$distances
    accepted <- all(!is.na(state$thresholds)) &&
      all(distances <= state$thresholds)

    stageStoredCollection(
      model_result$summaries, "summaries", state$store_summaries,
      accepted, attempt_id, 0L, spec$job_id,
      state$tmp_object_store_root
    )
    stageStoredCollection(
      model_result$outputs, "outputs", state$store_outputs,
      accepted, attempt_id, 0L, spec$job_id,
      state$tmp_object_store_root
    )

    particle_values <- proposed_particle[
      setdiff(names(proposed_particle), "job_id")
    ]
    rows[[attempt_index]] <- buildParticleRow(
      particle_values,
      list(
        attempt_id = attempt_id,
        job_id = spec$job_id,
        accepted = accepted,
        retained = FALSE,
        batch_id = spec$batch_id,
        attempt_index = attempt_index
      ),
      stats::setNames(as.list(distances), state$dist_names),
      c(state$column_names, "batch_id", "attempt_index")
    )
  }

  dplyr::bind_rows(rows)
}

runABCBatch <- function(kind = c("smc", "rejection"), spec, state_path,
                        batch_dir) {
  kind <- match.arg(kind)
  started_at <- Sys.time()
  control_path <- file.path(batch_dir, "result.rds")

  tryCatch({
    state <- new.env(parent = globalenv())
    load(state_path, envir = state)
    if (!is.null(state$model_def)) {
      sys.source(state$model_def, envir = state)
    }
    set.seed(spec$seed)
    particles <- if (kind == "smc") {
      runSMCBatch(spec, state)
    } else {
      runRejectionBatch(spec, state)
    }
    particle_path <- file.path(batch_dir, "particles.rds")
    writeRDSAtomically(particles, particle_path)
    result <- list(
      protocol_version = 1L,
      status = "ok",
      kind = kind,
      batch_id = spec$batch_id,
      generation = spec$generation,
      job_id = spec$job_id,
      n_attempted = nrow(particles),
      n_accepted = sum(particles$accepted),
      particles_path = normalizePath(
        particle_path, winslash = "/", mustWork = FALSE
      ),
      started_at = started_at,
      finished_at = Sys.time()
    )
    writeRDSAtomically(result, control_path)
    result
  }, error = function(error) {
    result <- list(
      protocol_version = 1L,
      status = "error",
      kind = kind,
      batch_id = spec$batch_id,
      generation = spec$generation,
      job_id = spec$job_id,
      n_attempted = 0L,
      n_accepted = 0L,
      message = conditionMessage(error),
      call = paste(deparse(conditionCall(error)), collapse = " "),
      started_at = started_at,
      finished_at = Sys.time()
    )
    try(writeRDSAtomically(result, control_path), silent = TRUE)
    result
  })
}

launchLocalBatch <- function(kind, spec, state_path, batch_root) {
  batch_dir <- file.path(
    batch_root,
    sprintf("generation_%04d", spec$generation),
    sprintf("batch_%08d", spec$batch_id)
  )
  dir.create(batch_dir, recursive = TRUE, showWarnings = FALSE)
  process <- callr::r_bg(
    func = runABCBatch,
    args = list(
      kind = kind,
      spec = spec,
      state_path = state_path,
      batch_dir = batch_dir
    ),
    stdout = file.path(batch_dir, "stdout.log"),
    stderr = file.path(batch_dir, "stderr.log"),
    supervise = TRUE,
    package = TRUE
  )
  list(process = process, spec = spec, batch_dir = batch_dir)
}

collectLocalBatch <- function(handle) {
  exit_status <- handle$process$get_exit_status()
  if (!identical(exit_status, 0L)) {
    stderr <- tryCatch(
      paste(handle$process$read_all_error_lines(), collapse = "\n"),
      error = function(...) ""
    )
    stop(sprintf(
      "Batch %d exited with status %s.%s",
      handle$spec$batch_id, exit_status,
      if (nzchar(stderr)) paste0("\n", stderr) else ""
    ), call. = FALSE)
  }
  result <- handle$process$get_result()
  if (!identical(result$status, "ok")) {
    stop(sprintf(
      "Batch %d failed: %s",
      handle$spec$batch_id, result$message %||% "unknown error"
    ), call. = FALSE)
  }
  result
}

runLocalBatches <- function(kind, generation, state_path, batch_root,
                            target_accepted, max_attempts, batch_size,
                            max_concurrent_jobs,
                            acceptance_rate_min = NULL,
                            lhs_indices = NULL, progress = NULL) {
  batch_size <- validateBatchSize(batch_size)
  max_concurrent_jobs <- validatePositiveInteger(
    max_concurrent_jobs, "max_concurrent_jobs"
  )
  max_attempts <- validatePositiveInteger(max_attempts, "max_attempts")
  state_path <- normalizePath(state_path, winslash = "/", mustWork = TRUE)
  batch_root <- normalizePath(batch_root, winslash = "/", mustWork = FALSE)

  active <- list()
  completed <- list()
  next_batch_id <- 1L
  attempts_reserved <- 0L
  attempts_completed <- 0L
  accepted_completed <- 0L
  stop_submitting <- FALSE
  stop_reason <- NULL

  terminate_active <- function() {
    for (handle in active) {
      if (handle$process$is_alive()) handle$process$kill()
    }
  }
  completed_normally <- FALSE
  on.exit(if (!completed_normally) terminate_active(), add = TRUE)

  repeat {
    while (!stop_submitting && length(active) < max_concurrent_jobs &&
           attempts_reserved < max_attempts) {
      n_attempts <- min(batch_size, max_attempts - attempts_reserved)
      batch_lhs_indices <- NULL
      if (!is.null(lhs_indices)) {
        remaining <- length(lhs_indices) - attempts_reserved
        if (remaining <= 0L) break
        n_attempts <- min(n_attempts, remaining)
        batch_lhs_indices <- lhs_indices[
          attempts_reserved + seq_len(n_attempts)
        ]
      }
      occupied_jobs <- vapply(
        active, function(x) x$spec$job_id, integer(1)
      )
      job_id <- setdiff(seq_len(max_concurrent_jobs), occupied_jobs)[[1L]]
      spec <- newBatchSpec(
        next_batch_id, generation, job_id, n_attempts, batch_lhs_indices
      )
      active[[as.character(next_batch_id)]] <- launchLocalBatch(
        kind, spec, state_path, batch_root
      )
      attempts_reserved <- attempts_reserved + n_attempts
      next_batch_id <- next_batch_id + 1L
    }

    if (!length(active)) break

    finished_ids <- names(active)[!vapply(
      active, function(x) x$process$is_alive(), logical(1)
    )]
    if (!length(finished_ids)) {
      Sys.sleep(0.1)
      next
    }

    for (id in finished_ids) {
      result <- collectLocalBatch(active[[id]])
      particles <- readRDS(result$particles_path)
      completed[[length(completed) + 1L]] <- particles
      attempts_completed <- attempts_completed + result$n_attempted
      accepted_completed <- accepted_completed + result$n_accepted
      active[[id]] <- NULL
    }

    if (is.function(progress)) {
      progress(attempts_completed, accepted_completed)
    }
    if (!is.null(target_accepted) &&
        accepted_completed >= target_accepted) {
      stop_submitting <- TRUE
      stop_reason <- "target_reached"
    }
    if (!is.null(acceptance_rate_min) && acceptance_rate_min > 0 &&
        attempts_completed > 0L &&
        attempts_completed >= ceiling(1 / acceptance_rate_min) &&
        accepted_completed / attempts_completed < acceptance_rate_min) {
      stop_submitting <- TRUE
      stop_reason <- "acceptance_rate"
    }
    if (attempts_reserved >= max_attempts && !length(active)) {
      stop_reason <- stop_reason %||% "max_attempts"
    }
  }

  completed_normally <- TRUE
  particles <- dplyr::bind_rows(completed)
  if (nrow(particles)) {
    particles <- particles[
      order(particles$batch_id, particles$attempt_index), , drop = FALSE
    ]
    rownames(particles) <- NULL
  }
  list(
    particles = particles,
    attempts = attempts_completed,
    accepted = accepted_completed,
    stop_reason = stop_reason %||% "completed"
  )
}

submitClusterWave <- function(cluster_type, template, first_task, last_task,
                              max_concurrent_jobs, stdout_dir, stderr_dir,
                              script_path) {
  cluster_script <- sprintf(
    template, first_task, last_task, max_concurrent_jobs,
    stdout_dir, stderr_dir, script_path
  )
  launcher_path <- file.path(dirname(script_path), "launch.sh")
  writeLines(cluster_script, launcher_path)
  command <- if (cluster_type == "slurm") "sbatch" else "qsub"
  output <- system2(command, launcher_path, stdout = TRUE, stderr = TRUE)
  status <- attr(output, "status") %||% 0L
  if (status != 0L) {
    stop(sprintf(
      "Unable to submit the %s batch wave:\n%s",
      cluster_type, paste(output, collapse = "\n")
    ), call. = FALSE)
  }
  matches <- regmatches(output, gregexpr("[0-9]+", output))
  identifiers <- unlist(matches, use.names = FALSE)
  if (!length(identifiers)) {
    stop("Unable to determine the submitted cluster job identifier.",
         call. = FALSE)
  }
  if (cluster_type == "slurm") {
    utils::tail(identifiers, 1L)
  } else {
    identifiers[[1L]]
  }
}

clusterJobIsActive <- function(cluster_type, job_id) {
  if (cluster_type == "slurm") {
    output <- suppressWarnings(system2(
      "squeue", c("-h", "-j", job_id), stdout = TRUE, stderr = FALSE
    ))
    length(output) > 0L
  } else {
    output <- suppressWarnings(system2(
      "qstat", c("-j", job_id), stdout = TRUE, stderr = FALSE
    ))
    identical(attr(output, "status") %||% 0L, 0L)
  }
}

runClusterBatches <- function(kind, generation, state_path, batch_root,
                              target_accepted, max_attempts, batch_size,
                              max_concurrent_jobs, cluster_type,
                              slurm_script_template, sge_script_template,
                              stdout_dir, stderr_dir,
                              acceptance_rate_min = NULL,
                              lhs_indices = NULL, progress = NULL) {
  batch_size <- validateBatchSize(batch_size)
  if (!cluster_type %in% c("slurm", "sge")) {
    stop("Cluster type not supported in the current version.", call. = FALSE)
  }
  max_concurrent_jobs <- validatePositiveInteger(
    max_concurrent_jobs, "max_concurrent_jobs"
  )
  max_attempts <- validatePositiveInteger(max_attempts, "max_attempts")
  state_path <- normalizePath(state_path, winslash = "/", mustWork = TRUE)
  batch_root <- normalizePath(batch_root, winslash = "/", mustWork = FALSE)
  dir.create(stdout_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(stderr_dir, recursive = TRUE, showWarnings = FALSE)

  completed <- list()
  next_batch_id <- 1L
  attempts_completed <- 0L
  accepted_completed <- 0L
  stop_reason <- NULL
  wave_id <- 1L

  repeat {
    if (attempts_completed >= max_attempts) {
      stop_reason <- "max_attempts"
      break
    }
    if (!is.null(target_accepted) && accepted_completed >= target_accepted) {
      stop_reason <- "target_reached"
      break
    }

    remaining_attempts <- max_attempts - attempts_completed
    if (!is.null(lhs_indices)) {
      remaining_attempts <- min(
        remaining_attempts, length(lhs_indices) - attempts_completed
      )
      if (remaining_attempts <= 0L) break
    }
    wave_jobs <- min(
      max_concurrent_jobs,
      ceiling(remaining_attempts / batch_size)
    )
    specs <- vector("list", wave_jobs)
    batch_dirs <- vector("list", wave_jobs)
    reserved_in_wave <- 0L
    for (task_id in seq_len(wave_jobs)) {
      n_attempts <- min(
        batch_size, remaining_attempts - reserved_in_wave
      )
      batch_lhs_indices <- NULL
      if (!is.null(lhs_indices)) {
        batch_lhs_indices <- lhs_indices[
          attempts_completed + reserved_in_wave + seq_len(n_attempts)
        ]
      }
      specs[[task_id]] <- newBatchSpec(
        next_batch_id, generation, task_id, n_attempts, batch_lhs_indices
      )
      batch_dirs[[task_id]] <- file.path(
        batch_root,
        sprintf("generation_%04d", generation),
        sprintf("batch_%08d", next_batch_id)
      )
      dir.create(batch_dirs[[task_id]], recursive = TRUE, showWarnings = FALSE)
      reserved_in_wave <- reserved_in_wave + n_attempts
      next_batch_id <- next_batch_id + 1L
    }

    wave_dir <- file.path(
      batch_root, sprintf("generation_%04d", generation),
      sprintf("wave_%06d", wave_id)
    )
    dir.create(wave_dir, recursive = TRUE, showWarnings = FALSE)
    manifest_path <- file.path(wave_dir, "manifest.rds")
    writeRDSAtomically(
      list(specs = specs, batch_dirs = batch_dirs), manifest_path
    )
    script_path <- file.path(wave_dir, "run_batch.R")
    script <- c(
      "library(BRREWABC)",
      "args <- commandArgs(trailingOnly = TRUE)",
      "task_id <- as.integer(args[[1L]])",
      sprintf("manifest <- readRDS(%s)", encodeString(
        normalizePath(manifest_path, winslash = "/", mustWork = TRUE),
        quote = "\""
      )),
      sprintf("state_path <- %s", encodeString(state_path, quote = "\"")),
      sprintf("kind <- %s", encodeString(kind, quote = "\"")),
      "result <- BRREWABC:::runABCBatch(kind, manifest$specs[[task_id]], state_path, manifest$batch_dirs[[task_id]])",
      "if (!identical(result$status, \"ok\")) quit(status = 1L)"
    )
    writeLines(script, script_path)
    template <- if (cluster_type == "slurm") {
      slurm_script_template
    } else {
      sge_script_template
    }
    job_id <- submitClusterWave(
      cluster_type, template, 1L, wave_jobs, max_concurrent_jobs,
      stdout_dir, stderr_dir, script_path
    )

    result_paths <- file.path(unlist(batch_dirs), "result.rds")
    job_seen_active <- FALSE
    inactive_checks <- 0L
    repeat {
      complete <- file.exists(result_paths)
      if (all(complete)) break
      Sys.sleep(0.25)
      active <- clusterJobIsActive(cluster_type, job_id)
      job_seen_active <- job_seen_active || active
      inactive_checks <- if (active) 0L else inactive_checks + 1L
      if ((job_seen_active && !active) || inactive_checks >= 20L) {
        Sys.sleep(0.5)
        if (!all(file.exists(result_paths))) {
          missing <- vapply(
            specs[!file.exists(result_paths)],
            function(x) x$batch_id, integer(1)
          )
          stop(sprintf(
            "Cluster job %s ended without results for batch(es): %s.",
            job_id, paste(missing, collapse = ", ")
          ), call. = FALSE)
        }
      }
    }

    for (result_path in result_paths) {
      result <- readRDS(result_path)
      if (!identical(result$status, "ok")) {
        stop(sprintf(
          "Batch %d failed: %s",
          result$batch_id, result$message %||% "unknown error"
        ), call. = FALSE)
      }
      particles <- readRDS(result$particles_path)
      completed[[length(completed) + 1L]] <- particles
      attempts_completed <- attempts_completed + result$n_attempted
      accepted_completed <- accepted_completed + result$n_accepted
    }
    if (is.function(progress)) {
      progress(attempts_completed, accepted_completed)
    }
    if (!is.null(acceptance_rate_min) && acceptance_rate_min > 0 &&
        attempts_completed > 0L &&
        attempts_completed >= ceiling(1 / acceptance_rate_min) &&
        accepted_completed / attempts_completed < acceptance_rate_min) {
      stop_reason <- "acceptance_rate"
      break
    }
    wave_id <- wave_id + 1L
  }

  particles <- dplyr::bind_rows(completed)
  if (nrow(particles)) {
    particles <- particles[
      order(particles$batch_id, particles$attempt_index), , drop = FALSE
    ]
    rownames(particles) <- NULL
  }
  list(
    particles = particles,
    attempts = attempts_completed,
    accepted = accepted_completed,
    stop_reason = stop_reason %||% "completed"
  )
}

#' Run a batch task of ABC-SMC
#'
#' This compatibility wrapper runs a batch specification stored in
#' `batch_manifest_path` inside the saved ABC state.
#' @param job_id array-task index in the batch manifest.
#' @param path_to_abc_state path to the saved ABC state.
#' @return the batch control object, invisibly.
#' @export
subjob_smc <- function(job_id, path_to_abc_state) {
  state <- new.env(parent = globalenv())
  load(path_to_abc_state, envir = state)
  if (is.null(state$batch_manifest_path)) {
    stop("The ABC state does not contain a batch manifest.", call. = FALSE)
  }
  manifest <- readRDS(state$batch_manifest_path)
  invisible(runABCBatch(
    "smc", manifest$specs[[job_id]], path_to_abc_state,
    manifest$batch_dirs[[job_id]]
  ))
}

#' Run a batch task of ABC rejection
#'
#' This compatibility wrapper runs a batch specification stored in
#' `batch_manifest_path` inside the saved ABC state.
#' @inheritParams subjob_smc
#' @return the batch control object, invisibly.
#' @export
subjob_rejection <- function(job_id, path_to_abc_state) {
  state <- new.env(parent = globalenv())
  load(path_to_abc_state, envir = state)
  if (is.null(state$batch_manifest_path)) {
    stop("The ABC state does not contain a batch manifest.", call. = FALSE)
  }
  manifest <- readRDS(state$batch_manifest_path)
  invisible(runABCBatch(
    "rejection", manifest$specs[[job_id]], path_to_abc_state,
    manifest$batch_dirs[[job_id]]
  ))
}
