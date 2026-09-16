
#' Run ABC rejection inference in parallel
#'
#' @param model_list a list linking model name ( character string) to
#' associated function
#' @param model_def a R file containing only the model(s) function(s)
#' @param prior_dist a list linking model name (character string) to a list
#' describing the prior distribution of each parameter to be estimated
#' @param ss_obs the observed summary statistics
#' @param nb_acc_prtcl the number of particles (per model, the total
#' number corresponding to this number multiplied by the number of models) to
#' be accepted
#' @param thresholds a value of the threshold to be used to select acceptable
#' particles (currently, using multiple distances is not yet supported for this
#' method, but will be considered in the future). If NA, no particle will be
#' accepted and a number of parricles equal to max_attempts will be generated
#' @param max_attempts the maximum number of particles to be tested during an
#' iteration, beyond which the procedure stops (in order to prevent
#' excessively long computations)
#' @param acceptance_rate_min the acceptance rate below which the procedure
#' stops (in order to prevent excessively long computations)
#' @param experiment_folderpath the folder in which to carry out the estimation
#' procedure and save the results
#' @param on_cluster whether or not the procedure is run on a
#' computation cluster
#' @param cluster_type cluster type used (sge and slurm currently supported)
#' @param slurm_script_template script used to launch jobs on a slurm cluster
#' @param sge_script_template script used to launch jobs on a sge cluster
#' @param max_concurrent_jobs maximum number of jobs/tasks run in parallel
#' @param batch_size number of simulations performed by a local or cluster
#' worker before returning control to the coordinator. Process startup and
#' package loading occur once per batch, so use a sufficiently large value for
#' inexpensive simulations.
#' @param store_summaries which summary statistics to store in Parquet files:
#' `"none"`, `"retained"`, `"accepted"`, or `"all"`
#' @param store_outputs which model outputs to store in Parquet files, using the
#' same policies as `store_summaries`
# #' @param abc_user_param_file_path an R file containing the algorithm's
# #' parameters (usage not recommended, included in this version for reasons of
# #' compatibility with the procedure in script form used in some projects)
#' @param verbose whether or not to display specific information
#' @param progressbar whether or not to display progressbar
#'
#' @return a list containing accepted particles, all tested particles, and a
#' `storage` descriptor for summary statistics and model outputs in Parquet.
#' @export
#' @include subjob.R saveEnvir.R
#'
#' @examples
#' library(BRREWABC)
abcrejection <- function(model_list = list(), # required
                         model_def = NULL,
                         prior_dist = list(), # required
                         ss_obs = NA, # required
                         nb_acc_prtcl = 1000,
                         thresholds = NA,
                         max_attempts = 100000,
                         acceptance_rate_min = 0.01,
                         experiment_folderpath = "./",
                         on_cluster = FALSE,
                         cluster_type = NULL, # "slurm" # "sge"
                         slurm_script_template = '#!/bin/bash
# THE FOLLOWING SECTION SHOULD NOT BE MODIFIED
#SBATCH --job-name=job-array_%%A_%%a   # nom du job
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --hint=nomultithread
#SBATCH --time=24:00:00
#SBATCH --array=%s-%s%%%d
output_fpath=%s
error_fpath=%s
#SBATCH --output=$output_fpath/output_%%A_%%a.out
#SBATCH --error=$error_fpath/error_%%A_%%a.out
mkdir -p $output_fpath
mkdir -p $error_fpath
Rscript %s $SLURM_ARRAY_TASK_ID
', # TODO : queue selection via a function argument
                         sge_script_template = '#!/bin/bash
#$ -S /bin/bash
#$ -N subjob_abcrejection_prlll
# #$ -q "short.q|long.q"
# THE FOLLOWING SECTION SHOULD NOT BE MODIFIED
#$ -cwd
#$ -V
#$ -t %s-%s
#$ -tc %d
#$ -o /dev/null
#$ -e /dev/null
output_fpath=%s
error_fpath=%s
mkdir -p $output_fpath
mkdir -p $error_fpath
Rscript %s $SGE_TASK_ID >$output_fpath/subjob.${SGE_TASK_ID}.out 2>$error_fpath/subjob.${SGE_TASK_ID}.err
', # TODO : queue selection via a function argument
                         max_concurrent_jobs = 1,
                         batch_size = 50,
                         store_summaries = c("none", "retained", "accepted", "all"),
                         store_outputs = c("none", "retained", "accepted", "all"),
                         verbose = FALSE,
                         progressbar = FALSE) {

  store_summaries <- normalizeStoragePolicy(store_summaries)
  store_outputs <- normalizeStoragePolicy(store_outputs)

  tmp_folder_path <- "tmp"
  results_folder_path <- "res"
  if (experiment_folderpath != "") {
    tmp_folder_path <- file.path(experiment_folderpath, tmp_folder_path)
    results_folder_path <- file.path(experiment_folderpath, results_folder_path)
  }

  tmp_local_task_std_out <- file.path(tmp_folder_path, "std_out")
  tmp_local_task_std_err <- file.path(tmp_folder_path, "std_err")
  tmp_current_abc_state <- file.path(tmp_folder_path, "currentABCState.RData")
  tmp_object_store_root <- file.path(tmp_folder_path, "stored_objects")
  tmp_batch_root <- file.path(tmp_folder_path, "batches")

  results_folder_path_CSV <- file.path(results_folder_path, "csv")
  results_folder_path_FIGS <- file.path(results_folder_path, "figs")
  storage_root <- file.path(results_folder_path, "parquet")
  accepted_particles_filepath <- file.path(results_folder_path_CSV, "all_accepted_particles.csv")
  all_tested_particles_filepath <- file.path(results_folder_path_CSV, "all_particles.csv")

  #

  for (folder_path in c(tmp_folder_path, results_folder_path, results_folder_path_CSV, results_folder_path_FIGS, storage_root)) {
    if (verbose) {cat(paste0("Check folder_path for : ", folder_path, "\n"))}
    # Check if the folder path exists
    if (!dir.exists(folder_path)) {
      # Folder does not exist, create the folder
      if (dir.create(folder_path, recursive = TRUE)) {
        if (verbose) {cat("Folder created successfully.\n")}
      } else {
        if (verbose) {cat("Error: Unable to create folder.\n")}
      }
    } else {
      if (verbose) {cat("Folder already exists.\n")}
    }
  }
  unlink(tmp_object_store_root, recursive = TRUE)
  dir.create(tmp_object_store_root, recursive = TRUE, showWarnings = FALSE)

  nb_threshold <- length(thresholds)
  dist_names <- paste0("dist", as.character(seq(1, nb_threshold, 1)))
  model_names <- names(model_list)
  param_names <- unique(Reduce(c, sapply(prior_dist, function(x) sapply(x, `[[`, 1))))
  column_names <- c("attempt_id", "job_id", "accepted", "retained",
                    "model", param_names, dist_names)

  # define the total number of particles to accept before next gen, that will be
  # used as a upper limit for the number of simulation to run (avoid a while
  # loop without stopping criterion)
  tot_nb_acc_prtcl <- nb_acc_prtcl * length(model_names)

  acc_particles <- stats::setNames(data.frame(matrix(ncol = length(column_names), nrow = 0), stringsAsFactors=FALSE), column_names)
  all_tested_particles <- stats::setNames(data.frame(matrix(ncol = length(column_names), nrow = 0), stringsAsFactors=FALSE), column_names)

  utils::write.csv(acc_particles, accepted_particles_filepath, row.names=FALSE, quote=FALSE)
  utils::write.csv(all_tested_particles, all_tested_particles_filepath, row.names=FALSE, quote=FALSE)

  # print(ls()) # DEBUG
  unlink(file.path(storage_root, "manifest.parquet"))
  var_to_save <- c( "model_def", "model_list", "prior_dist", "ss_obs", "max_concurrent_jobs", "accepted_particles_filepath", "all_tested_particles_filepath", "tmp_object_store_root", "store_summaries", "store_outputs", "dist_names", "model_names", "param_names", "column_names", "tot_nb_acc_prtcl", "thresholds")
    # saveEnvir(var_to_save, tmp_current_abc_state) # TODO : not working, need to fix this
  do.call("save", c(var_to_save, list(file = tmp_current_abc_state)))

  pb <- NULL
  if (progressbar) {
    pb <- progress::progress_bar$new(
      format = "[:bar] :percent (ar: :accrate | :nbattempt) | eta: :eta (:elapsed)",
      clear = FALSE,
      total = tot_nb_acc_prtcl
    )
  }
  update_progress <- function(attempted, accepted) {
    if (!is.null(pb)) {
      rate <- if (attempted) accepted / attempted else 0
      pb$update(
        min(accepted / tot_nb_acc_prtcl, 1),
        tokens = list(
          accrate = format(round(rate, 3), nsmall = 3),
          nbattempt = attempted
        )
      )
    }
  }
  target_accepted <- if (all(is.na(thresholds))) NULL else tot_nb_acc_prtcl
  acceptance_limit <- if (all(is.na(thresholds))) NULL else acceptance_rate_min
  elapsed <- system.time({
    batch_result <- if (on_cluster) {
      runClusterBatches(
        "rejection", 0L, tmp_current_abc_state, tmp_batch_root,
        target_accepted, max_attempts, batch_size, max_concurrent_jobs,
        cluster_type, slurm_script_template, sge_script_template,
        tmp_local_task_std_out, tmp_local_task_std_err,
        acceptance_limit, progress = update_progress
      )
    } else {
      runLocalBatches(
        "rejection", 0L, tmp_current_abc_state, tmp_batch_root,
        target_accepted, max_attempts, batch_size, max_concurrent_jobs,
        acceptance_limit, progress = update_progress
      )
    }
  })
  all_tested_particles <- batch_result$particles
  accepted_indices <- which(all_tested_particles$accepted)
  if (!is.null(target_accepted)) {
    accepted_indices <- utils::head(accepted_indices, target_accepted)
  }
  all_tested_particles$retained <- FALSE
  all_tested_particles$retained[accepted_indices] <- TRUE
  acc_particles <- all_tested_particles[accepted_indices, , drop = FALSE]
  internal_columns <- c("batch_id", "attempt_index")
  all_tested_particles <- all_tested_particles[
    setdiff(names(all_tested_particles), internal_columns)
  ]
  acc_particles <- acc_particles[setdiff(names(acc_particles), internal_columns)]
  utils::write.csv(acc_particles, accepted_particles_filepath,
                   row.names = FALSE, quote = FALSE)
  utils::write.csv(all_tested_particles, all_tested_particles_filepath,
                   row.names = FALSE, quote = FALSE)
  if (verbose) {
    cat(sprintf(
      "Computation time - user : %.3f s | system : %.3f s | elapsed : %.3f s \n",
      elapsed["user.self"],
      elapsed["sys.self"],
      elapsed["elapsed"]
    ))
  }

  #
  retained_ids <- acc_particles$attempt_id
  persistStoredGeneration(tmp_object_store_root, storage_root, 0L,
                          all_tested_particles$attempt_id, retained_ids,
                          store_summaries, store_outputs)

  if (verbose) {
    cat("Experiment done!", "\n")
  }
  # cleaning
  unlink(tmp_folder_path, recursive = TRUE)
  #
  storage <- list(
    path = normalizePath(storage_root, winslash = "/", mustWork = FALSE),
    format = "parquet",
    manifest = file.path(normalizePath(storage_root, winslash = "/", mustWork = FALSE), "manifest.parquet")
  )
  return(list("acc_particles" = acc_particles,
              "all_tested_particles" = all_tested_particles,
              "storage" = storage))
}
