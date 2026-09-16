
#' Run ABC-SMC inference in parallel
#'
#' @param model_list a list linking model name ( character string) to
#' associated function
#' @param model_def a R file containing only the model(s) function(s)
#' @param prior_dist a list linking model name (character string) to a list
#' describing the prior distribution of each parameter to be estimated
#' @param ss_obs the observed summary statistics
#' @param max_number_of_gen the maximum number of generations to be performed
#' @param nb_acc_prtcl_per_gen the number of particles (per model, the total
#' number corresponding to this number multiplied by the number of models) to
#' be accepted during a generation before moving on to the next one
#' @param var_prtrbtn_krnl the standard deviation to be used for the
#' perturbation kernel, if NA the empirical standard deviation will
#' be used instead
#' @param var_prtrbtn_krnl_min the minimum standard deviation to be used for
#' the perturbation kernel (must be greater than zero)
#' @param scaling_sd_prtrbtn_krnl scaling parameter for increasing or decreasing
#' the value of the standard deviation used for the perturbation kernel
#' @param model_jump_prob probability of changing model when creating
#' a new particle
#' @param nb_threshold number of thresholds used, corresponding to the number
#' of distances returned by the function(s) stipulated in the
#' model list (model_list)
#' @param new_threshold_quantile the quantile to be used to decrease the
#' threshold(s) during iterations
#' @param distance_threshold_min the minimum threshold below which the procedure
#' stops (in order to prevent excessively long computations)
#' @param max_attempts the maximum number of particles to be tested during an
#' iteration, beyond which the procedure stops (in order to prevent
#' excessively long computations)
#' @param acceptance_rate_min the acceptance rate below which the procedure
#' stops (in order to prevent excessively long computations)
#' @param use_lhs_for_first_iter whether or not to use a LHS sampling method
#' for the first iteration
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
#' @param storage_chunk_rows maximum number of stored summary/output rows kept
#' in a worker buffer before an atomic Parquet fragment is written.
#' @param storage_chunk_mb approximate maximum buffer size, in MiB, before an
#' atomic Parquet fragment is written.
# #' @param abc_user_param_file_path an R file containing the algorithm's
# #' parameters (usage not recommended, included in this version for reasons of
# #' compatibility with the procedure in script form used in some projects)
#' @param previous_gens an object (dataframe) containing previous results (set
#' of iterations), in order to start from the last iteration performed
#' @param previous_epsilons an object (dataframe) containing previous results
#' (set of thresholds), in order to start from the last iteration performed
#' @param store_summaries which summary statistics to store in Parquet files:
#' `"none"`, `"retained"`, `"accepted"`, or `"all"`
#' @param store_outputs which model outputs to store in Parquet files, using the
#' same policies as `store_summaries`
#' @param verbose whether or not to display specific information
#' @param progressbar whether or not to display progressbar
#'
#' @return a list containing accepted particles, thresholds, and a `storage`
#' descriptor for summary statistics and model outputs stored in Parquet.
#' @export
#' @include createLHSfromPrior.R defineNextThreshold.R setEmpiricalSD.R subjob.R saveEnvir.R
#'
#' @examplesIf interactive()
#' library(BRREWABC)
#'
#' tmp_dir <- tempdir()
#'
#' # model definition
#' compute_dist = function(x, ss_obs){
#'     ss_sim = c( x[["alpha"]] + x[["beta"]] + rnorm(1,0,0.1),
#'            x[["alpha"]] * x[["beta"]] + rnorm(1,0,0.1) ) # a simple toy model
#'     dist = sum((ss_sim-ss_obs)^2)
#'     return(c(dist))
#' }
#'
#' MODEL_LIST <- list("m1" = compute_dist)
#' PRIOR_DIST <- list("m1" = list(c('alpha', 'unif', 0, 4),
#'                                c('beta', 'unif', 0, 1)))
#'
#' # create a reference trajectory
#' sum_stat_obs = c(2.0,0.75)
#'
#' # run abc smc procedure
#' res = abcsmc(model_list = MODEL_LIST, prior_dist = PRIOR_DIST,
#' ss_obs = sum_stat_obs, max_number_of_gen = 20, nb_acc_prtcl_per_gen = 2000,
#' new_threshold_quantile = 0.8, experiment_folderpath = tmp_dir,
#' max_concurrent_jobs = 2, verbose = FALSE)
#'
#' # get results and plots
#' all_accepted_particles = res$particles
#' all_thresholds = res$thresholds
#' plot_abcsmc_res(data = all_accepted_particles, prior = PRIOR_DIST, colorpal = "YlOrBr", filename = file.path(tmp_dir, "abcsmc_results.png"))
#' plot_thresholds(data = all_thresholds, nb_threshold = 1, colorpal = "YlOrBr", filename = file.path(tmp_dir, "thresholds.png"))
#' plot_ess(data = all_accepted_particles, colorpal = "YlOrBr", filename = file.path(tmp_dir, "ess.png"))
#' plot_densityridges(data = all_accepted_particles, prior = PRIOR_DIST, colorpal = "YlOrBr", filename = file.path(tmp_dir, "densityridges.png"))
abcsmc <- function(model_list = list(), # required
                   model_def = NULL,
                   prior_dist = list(), # required
                   ss_obs = NA, # required
                   max_number_of_gen = 15,
                   nb_acc_prtcl_per_gen = 1000,
                   var_prtrbtn_krnl = NA,
                   var_prtrbtn_krnl_min = 0.01,
                   scaling_sd_prtrbtn_krnl = 1.0,
                   model_jump_prob = 0.1,
                   nb_threshold = 1,
                   new_threshold_quantile = 0.8,
                   distance_threshold_min = 0.01,
                   max_attempts = 100000,
                   acceptance_rate_min = 0.01,
                   use_lhs_for_first_iter = TRUE,
                   experiment_folderpath = "./",
                   on_cluster = FALSE,
                   cluster_type = NULL, # "slurm" # "sge"
                   slurm_script_template = '#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --hint=nomultithread
#SBATCH --time=24:00:00
# THE FOLLOWING SECTION SHOULD NOT BE MODIFIED
#SBATCH --job-name=job-array_%%A_%%a   # nom du job
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
#$ -N subjob_abcsmc_prlll
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
                   storage_chunk_rows = 1000000,
                   storage_chunk_mb = 128,
                   # abc_user_param_file_path = NULL,
                   previous_gens = NA,
                   previous_epsilons = NA,
                   store_summaries = c("none", "retained", "accepted", "all"),
                   store_outputs = c("none", "retained", "accepted", "all"),
                   verbose = FALSE,
                   progressbar = FALSE) {

  store_summaries <- normalizeStoragePolicy(store_summaries)
  store_outputs <- normalizeStoragePolicy(store_outputs)

  # # first, load user param file if exist
  # if (!is.null(abc_user_param_file_path)) {
  #   source(abc_user_param_file_path)
  # }

  #

  tmp_folder_path <- "tmp"
  results_folder_path <- "res"
  if (experiment_folderpath != "") {
    tmp_folder_path <- file.path(experiment_folderpath, tmp_folder_path)
    results_folder_path <- file.path(experiment_folderpath, results_folder_path)
  }

  tmp_object_store_root <- file.path(tmp_folder_path, "stored_objects")
  tmp_local_task_std_out <- file.path(tmp_folder_path, "std_out")
  tmp_local_task_std_err <- file.path(tmp_folder_path, "std_err")
  tmp_current_abc_state <- file.path(tmp_folder_path, "currentABCState.RData")
  tmp_batch_root <- file.path(tmp_folder_path, "batches")

  results_folder_path_CSV <- file.path(results_folder_path, "csv")
  results_folder_path_FIGS <- file.path(results_folder_path, "figs")
  storage_root <- file.path(results_folder_path, "parquet")
  accepted_particles_filepath <- file.path(results_folder_path_CSV, "all_accepted_particles.csv")
  last_accepted_particles_filepath <- file.path(results_folder_path_CSV, "last_accepted_particles.csv")
  all_particles_filepath <- file.path(results_folder_path_CSV, "all_particles.csv")
  thresholds_filepath <- file.path(results_folder_path_CSV, "thresholds.csv")
  lhs_used_for_first_gen_filepath <- file.path(results_folder_path_CSV, "lhs_for_first_gen.csv")

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

  #


  dist_names <- paste0("dist", as.character(seq(1, nb_threshold, 1)))
  model_names <- names(model_list)
  param_names <- unique(Reduce(c, sapply(prior_dist, function(x) sapply(x, `[[`, 1))))
  column_names <- c("gen", "attempt_id", "job_id", "accepted", "retained",
                    "model", param_names, "pWeight", dist_names)

  # define the total number of particles to accept before next gen, that will be
  # used as a upper limit for the number of simulation to run (avoid a while
  # loop without stopping criterion)
  nb_acc_prtcl_before_next_gen <- nb_acc_prtcl_per_gen * length(model_names)

  # define the vector of thresholds
  epsilon <- c()

  # define the vector of sd to use in the perturbation kernel
  empirical_sd <- list()
  for (mm in model_names) {
    empirical_sd[[mm]] <- list()
    for (pp in prior_dist[[mm]]) {
      empirical_sd[[mm]][[pp[1]]] <- var_prtrbtn_krnl
    }
  }

  if (length(scaling_sd_prtrbtn_krnl) == 1) {
    scaling_sd_prtrbtn_krnl_uniqval <- scaling_sd_prtrbtn_krnl
    scaling_sd_prtrbtn_krnl <- list()
    for (pp in param_names) {
      scaling_sd_prtrbtn_krnl[[pp]] <- scaling_sd_prtrbtn_krnl_uniqval
    }
  }

  all_prior_are_unif <- all(unique(Reduce(c, sapply(prior_dist, function(x) sapply(x, `[[`, 2)))) == "unif")

  #
  gen <- 1
  epsilon_not_improved <- 0

  all_acc_particles <- stats::setNames(data.frame(matrix(ncol = length(column_names), nrow = 0), stringsAsFactors=FALSE), column_names)
  previous_acc_particles <- stats::setNames(data.frame(matrix(ncol = length(column_names), nrow = 0), stringsAsFactors=FALSE), column_names)
  current_acc_particles <- stats::setNames(data.frame(matrix(ncol = length(column_names), nrow = 0), stringsAsFactors=FALSE), column_names)

  epsilons <- stats::setNames(data.frame(matrix(ncol = length(dist_names)+2, nrow = 0), stringsAsFactors=FALSE), c("gen", dist_names, "acceptance_rate"))

  lhs_first_gen = c()
  if (use_lhs_for_first_iter) {
    if (!all_prior_are_unif) {
      use_lhs_for_first_iter <- FALSE
    } else {
      lhs_first_gen = createLHSfromPrior(model_names, prior_dist, nb_acc_prtcl_per_gen)
      utils::write.csv(lhs_first_gen, lhs_used_for_first_gen_filepath, row.names=FALSE, quote=FALSE)
    }
  }

  # init with previous results if provided
  if (!is.null(dim(previous_gens))) {
    # all_acc_particles
    all_acc_particles <- previous_gens
    epsilons <- previous_epsilons
    gen <- max(all_acc_particles$gen)
    # previous_acc_particles : set based on last in all_acc_particles
    previous_acc_particles <- as.data.frame(all_acc_particles[all_acc_particles$gen == gen,])
    # current_acc_particles  # nothing to do
    # epsilon
    epsilon <- defineNextThreshold(previous_acc_particles, nb_threshold, dist_names, new_threshold_quantile)
    # calculate the empirical sd to be used in the perturbation kernel
    empirical_sd <- setEmpiricalSd(previous_acc_particles, model_names, prior_dist, var_prtrbtn_krnl, scaling_sd_prtrbtn_krnl, var_prtrbtn_krnl_min)
    # set the generation number
    gen <- gen + 1
  }
  if (is.null(dim(previous_gens))) unlink(file.path(storage_root, "manifest.parquet"))
  #
  while (gen <= max_number_of_gen) {
    nb_accepted <- 0
    totattempts <- 0
    if (verbose) {
      print(Sys.time())
      cat("gen", gen, "\n")
      cat("threshold:", epsilon, "\n")
      cat("prtrbtn_krnl_sd:", unlist(empirical_sd), "\n")
    }

    # print(ls()) # DEBUG
    var_to_save <- c( "model_def", "model_list", "prior_dist", "ss_obs", "model_jump_prob", "use_lhs_for_first_iter", "max_concurrent_jobs", "tmp_object_store_root", "store_summaries", "store_outputs", "storage_chunk_rows", "storage_chunk_mb", "dist_names", "model_names", "param_names", "column_names", "nb_acc_prtcl_before_next_gen", "epsilon", "empirical_sd", "lhs_first_gen", "gen", "previous_acc_particles")
      # saveEnvir(var_to_save, tmp_current_abc_state) # TODO : not working, need to fix this
    do.call("save", c(var_to_save, list(file = tmp_current_abc_state)))

    pb <- NULL
    if (progressbar) {
      pb <- progress::progress_bar$new(
        format = "gen :gen [:bar] :percent (ar: :accrate | :nbattempt) | eta: :eta (:elapsed)",
        clear = FALSE,
        total = nb_acc_prtcl_before_next_gen
      )
    }
    update_progress <- function(attempted, accepted) {
      if (!is.null(pb)) {
        rate <- if (attempted) accepted / attempted else 0
        updateProgressBar(
          pb,
          accepted / nb_acc_prtcl_before_next_gen,
          tokens = list(
            gen = gen,
            accrate = format(round(rate, 3), nsmall = 3),
            nbattempt = attempted
          )
        )
      }
    }
    lhs_indices <- if (gen == 1L && use_lhs_for_first_iter) {
      seq_len(nb_acc_prtcl_before_next_gen)
    } else {
      NULL
    }
    generation_max_attempts <- if (!is.null(lhs_indices)) {
      nb_acc_prtcl_before_next_gen
    } else {
      max_attempts
    }
    elapsed <- system.time({
      batch_result <- if (on_cluster) {
        runClusterBatches(
          "smc", gen, tmp_current_abc_state, tmp_batch_root,
          nb_acc_prtcl_before_next_gen, generation_max_attempts,
          batch_size, max_concurrent_jobs, cluster_type,
          slurm_script_template, sge_script_template,
          tmp_local_task_std_out, tmp_local_task_std_err,
          acceptance_rate_min, lhs_indices, update_progress
        )
      } else {
        runLocalBatches(
          "smc", gen, tmp_current_abc_state, tmp_batch_root,
          nb_acc_prtcl_before_next_gen, generation_max_attempts,
          batch_size, max_concurrent_jobs, acceptance_rate_min,
          lhs_indices, update_progress
        )
      }
    })
    tested_this_gen <- batch_result$particles
    totattempts <- batch_result$attempts
    nb_accepted <- batch_result$accepted
    current_acc_rate <- if (totattempts) nb_accepted / totattempts else 0
    current_iter_broke <- nb_accepted < nb_acc_prtcl_before_next_gen
    if (verbose) {
      print(Sys.time())
      cat("totattempts:", totattempts, "\n")
      cat("current_acc_rate:", format(round(current_acc_rate,digits=3),nsmall=3), "\n")
      cat(sprintf(
        "Computation time - user : %.3f s | system : %.3f s | elapsed : %.3f s \n",
        elapsed["user.self"],
        elapsed["sys.self"],
        elapsed["elapsed"]
      ))
    }
    retained_ids <- character()
    if (current_iter_broke == TRUE) {
      persistStoredGeneration(tmp_object_store_root, storage_root, gen,
                              tested_this_gen$attempt_id, retained_ids,
                              store_summaries, store_outputs,
                              batch_result$stored_fragments)
      break
    }
    #
    accepted_indices <- utils::head(
      which(tested_this_gen$accepted), nb_acc_prtcl_before_next_gen
    )
    tested_this_gen$retained <- FALSE
    tested_this_gen$retained[accepted_indices] <- TRUE
    current_acc_particles <- tested_this_gen[accepted_indices, , drop = FALSE]
    retained_ids <- current_acc_particles$attempt_id
    persistStoredGeneration(tmp_object_store_root, storage_root, gen,
                            tested_this_gen$attempt_id, retained_ids,
                            store_summaries, store_outputs,
                            batch_result$stored_fragments)
    internal_columns <- c("batch_id", "attempt_index")
    current_acc_particles <- current_acc_particles[
      setdiff(names(current_acc_particles), internal_columns)
    ]
    # normalise the weights
    for (mm in model_names) {
      if (mm %in% unique(current_acc_particles$model)) {
        current_acc_particles[current_acc_particles$model == mm,"pWeight"]= current_acc_particles[current_acc_particles$model == mm,"pWeight"]/ sum(current_acc_particles[current_acc_particles$model == mm,"pWeight"])
      }
    }
    utils::write.csv(current_acc_particles, last_accepted_particles_filepath, row.names=FALSE, quote=FALSE)
    # switch current and previous particles
    previous_acc_particles <- current_acc_particles
    # save current particles
    all_acc_particles = rbind(all_acc_particles, current_acc_particles)
    # clear current particles table for next generation
    current_acc_particles <- current_acc_particles[c(),]
    #
    if (gen > 1) { # save current epsilon if generation completed
      epsilons[nrow(epsilons) + 1,] <- c(gen, epsilon, current_acc_rate)
    }
    # save results
    utils::write.csv(all_acc_particles, accepted_particles_filepath, row.names=FALSE, quote=FALSE)
    utils::write.csv(epsilons, thresholds_filepath, row.names=FALSE, quote=FALSE)
    # STOPPING RULES
    # Stop the algo if there is no significant difference between two successive posterior distributions
    # if (gen > 1) {
    #     # To be implemented
    # }
    # Stop the algo if epsilon falls below the predetermined min threshold
    if ((gen > 1) && (all(epsilon <= distance_threshold_min))) {
      message("The distance threshold(s) (epsilon(s)) fall(s) below the predetermined min value!")
      print(epsilon)
      break
    }
    # calculate the empirical sd to be used in the perturbation kernel
    empirical_sd <- setEmpiricalSd(previous_acc_particles, model_names, prior_dist, var_prtrbtn_krnl, scaling_sd_prtrbtn_krnl, var_prtrbtn_krnl_min)
    # define the next threshold
    next_epsilon <- defineNextThreshold(previous_acc_particles, nb_threshold, dist_names, new_threshold_quantile)
    if ((gen > 1) && (all(next_epsilon == epsilon))) {
      epsilon_not_improved = epsilon_not_improved + 1
      if (epsilon_not_improved >= 2) {
        message("The distance threshold(s) (epsilon(s)) has not been improved over the last two iterations!")
        print(epsilons)
        break
      }
    }
    epsilon <- next_epsilon
    gen <- gen + 1
    if (verbose) {
      print(Sys.time())
      cat("-\n")
    }
  }
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
  return(list("particles" = all_acc_particles, "thresholds" = epsilons,
              "storage" = storage))
}
