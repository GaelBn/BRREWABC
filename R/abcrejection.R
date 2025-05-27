
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
# #' @param abc_user_param_file_path an R file containing the algorithm's
# #' parameters (usage not recommended, included in this version for reasons of
# #' compatibility with the procedure in script form used in some projects)
#' @param verbose whether or not to display specific information
#' @param progressbar whether or not to display progressbar
#'
#' @return a list containing two dataframes corresponding to (1) the particles
#' accepted and (2) all tested particles
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
                         verbose = FALSE,
                         progressbar = FALSE) {

  tmp_folder_path <- "tmp"
  results_folder_path <- "res"
  if (experiment_folderpath != "") {
    tmp_folder_path <- file.path(experiment_folderpath, tmp_folder_path)
    results_folder_path <- file.path(experiment_folderpath, results_folder_path)
  }

  tmp_local_task_std_out <- file.path(tmp_folder_path, "std_out")
  tmp_local_task_std_err <- file.path(tmp_folder_path, "std_err")
  tmp_current_abc_state <- file.path(tmp_folder_path, "currentABCState.RData")

  results_folder_path_CSV <- file.path(results_folder_path, "csv")
  results_folder_path_FIGS <- file.path(results_folder_path, "figs")
  accepted_particles_filepath <- file.path(results_folder_path_CSV, "all_accepted_particles.csv")
  all_tested_particles_filepath <- file.path(results_folder_path_CSV, "all_particles.csv")

  subjob_script_name <- "brrewabc_subjob.R"
  subjob_script_path <- file.path(tmp_folder_path, subjob_script_name)

  abc_rejection_array_job_script_path <- file.path(tmp_folder_path, "abc_rejection_array_job.sh")

  #

  for (folder_path in c(tmp_folder_path, results_folder_path, results_folder_path_CSV, results_folder_path_FIGS)) {
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

  #

  subjob_script_content <- "library(BRREWABC)\nargs <- commandArgs(trailingOnly = TRUE)\nid <- as.integer(args[1])\nsubjob_rejection(job_id = id, path_to_abc_state = '%s')\n"
  subjob_script <- sprintf(subjob_script_content, tmp_current_abc_state)
  writeLines(subjob_script, subjob_script_path)

  nb_threshold <- length(thresholds)
  dist_names <- paste0("dist", as.character(seq(1, nb_threshold, 1)))
  model_names <- names(model_list)
  param_names <- unique(Reduce(c, sapply(prior_dist, function(x) sapply(x, `[[`, 1))))
  column_names <- c("model", param_names, dist_names)

  # define the total number of particles to accept before next gen, that will be
  # used as a upper limit for the number of simulation to run (avoid a while
  # loop without stopping criterion)
  tot_nb_acc_prtcl <- nb_acc_prtcl * length(model_names)

  acc_particles <- stats::setNames(data.frame(matrix(ncol = length(column_names), nrow = 0), stringsAsFactors=FALSE), column_names)
  all_tested_particles <- stats::setNames(data.frame(matrix(ncol = length(column_names), nrow = 0), stringsAsFactors=FALSE), column_names)

  utils::write.csv(acc_particles, accepted_particles_filepath, row.names=FALSE, quote=FALSE)
  utils::write.csv(all_tested_particles, all_tested_particles_filepath, row.names=FALSE, quote=FALSE)

  #
  nb_accepted <- 0
  totattempts <- 0

  # print(ls()) # DEBUG
  var_to_save <- c( "model_def", "model_list", "prior_dist", "ss_obs", "max_concurrent_jobs", "accepted_particles_filepath", "all_tested_particles_filepath", "dist_names", "model_names", "param_names", "column_names", "tot_nb_acc_prtcl", "thresholds")
    # saveEnvir(var_to_save, tmp_current_abc_state) # TODO : not working, need to fix this
  do.call("save", c(var_to_save, list(file = tmp_current_abc_state)))

  #
  cluster_job_id <- NA
  local_processes <- c()
  if (on_cluster) {
    # submit a job array with as many job as tot_nb_acc_prtcl (for gen == 1, need at one job per each particle)
    cluster_script <- ""
    # Replace placeholders in the template with actual parameter values
    nbjobs <- max_concurrent_jobs
    if (cluster_type == "slurm") {
      cluster_script <- sprintf(slurm_script_template, 1, nbjobs, max_concurrent_jobs, tmp_local_task_std_out, tmp_local_task_std_err, subjob_script_path)
    } else if (cluster_type == "sge") {
      cluster_script <- sprintf(sge_script_template, 1, nbjobs, max_concurrent_jobs, tmp_local_task_std_out, tmp_local_task_std_err, subjob_script_path)
    } else {
      stop("Cluster type not supported in the current version")
    }
    # Write the modified cluster job script to a file
    writeLines(cluster_script, abc_rejection_array_job_script_path)
    # Submit the array job and capture the job ID
    if (cluster_type == "slurm") {
      cluster_job_id <- system(paste0('sbatch ', abc_rejection_array_job_script_path, ' | awk \'{print $4}\''), intern = TRUE)
      print(paste("Submitted Slurm array job with ID:", cluster_job_id))
    } else if (cluster_type == "sge") {
      cluster_job_id <- system(paste0('qsub ', abc_rejection_array_job_script_path, ' | awk -F "." \'{print $1}\' | awk \'{print $3}\''), intern = TRUE)
      print(paste("Submitted SGE array job with ID:", cluster_job_id))
    } else {
      stop("Cluster type not supported in the current version")
    }
  } else {
    # Parallel task on local machine
    for (task_index in 1:max_concurrent_jobs) {
      # Launch external script
      local_process <- callr::r_bg(
        func = function() {
          library(BRREWABC)
          subjob_rejection(job_id = task_index, path_to_abc_state = tmp_current_abc_state)
        },
        stdout = paste0(tmp_local_task_std_out,"abc_rejection_task_",task_index,".out"),
        stderr = paste0(tmp_local_task_std_err,"abc_rejection_task_",task_index,".err"),
        package = TRUE
      )
      local_processes <- c(local_processes, local_process)
    }
  }
  #
  # Create a progress bar
  if (progressbar && !on_cluster) {
    pb <- progress::progress_bar$new(
      format = "[:bar] :percent (ar: :accrate | :nbattempt) | eta: :eta (:elapsed)",
      clear = FALSE,
      total = tot_nb_acc_prtcl
    )
  }
  #
  elapsed <- system.time({
    while (nb_accepted < tot_nb_acc_prtcl) { # repeat until N particles accepted
      accepted_particles <- utils::read.csv(accepted_particles_filepath)
      all_tested_particles <- utils::read.csv(all_tested_particles_filepath)
      nb_accepted <- nrow(accepted_particles)
      totattempts <- nrow(all_tested_particles)
      current_acc_rate <- nb_accepted/totattempts
      # Increment the progress bar
      if (progressbar && !on_cluster) {
        pb$update(
          min(nb_accepted/tot_nb_acc_prtcl, 1.00),
          tokens = list(
            accrate = format(round(current_acc_rate,digits=3),nsmall=3),
            nbattempt = totattempts
          )
        )
      }
      if (totattempts > max_attempts) {
        message('\n', "The maximum number of attempts to accept the N particles has been reached!", '\n')
        break
      }
      if (all(!is.na(thresholds))) {
        if ((totattempts >= (1/acceptance_rate_min)) && (current_acc_rate < acceptance_rate_min)) {
          message('\n', "The acceptance rate has become too low (< specified acceptance_rate_min), the algorithm stops!")
          break
        }
      }
    }
  })
  # cancel subjob
  if (on_cluster) {
    # Cancel the array job
    if (cluster_type == "slurm") {
      cancel_command <- paste("scancel", cluster_job_id)
      system(cancel_command, intern = TRUE)
    } else if (cluster_type == "sge") {
      cancel_command <- paste("qdel", cluster_job_id)
      system(cancel_command, intern = TRUE)
    } else {
      stop("Cluster type not supported in the current version")
    }
  } else {
    # kill local process if needed
    for (local_process in local_processes) {
      if (local_process$is_alive()) {
        local_process$kill()
      }
    }
  }
  # Close the progress bar
  # pb$terminate()
  # cat('\n')
  #
  if (verbose) {
    message(sprintf(
      "Computation time - user : %.3f s | system : %.3f s | elapsed : %.3f s",
      elapsed["user.self"],
      elapsed["sys.self"],
      elapsed["elapsed"]
    ))
  }
  #
  acc_particles <- utils::read.csv(accepted_particles_filepath)
  all_tested_particles <- utils::read.csv(all_tested_particles_filepath)

  # acc_particles <- acc_particles[1:min(nrow(acc_particles), tot_nb_acc_prtcl),] # keep only the number of particle needed # TODO : improve comment
  # utils::write.csv(acc_particles, accepted_particles_filepath, row.names=FALSE, quote=FALSE)

  if (verbose) {
    cat("Experiment done!", "\n")
  }
  # cleaning
  if (on_cluster) {
    unlink(subjob_script_path)
    unlink(abc_rejection_array_job_script_path)
  }
  unlink(tmp_folder_path, recursive = TRUE)
  unlink(paste0(all_tested_particles_filepath, ".lck"))
  unlink(paste0(accepted_particles_filepath, ".lck"))
  #
  return(list("acc_particles" = acc_particles, "all_tested_particles" = all_tested_particles))
}
