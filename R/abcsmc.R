
#' Run ABC-SMC inference in parallel
#'
#' @param model_def a R file containing only the model(s) function(s)
#' @param model_list a list linking model name ( character string) to associated function
#' @param prior_dist a list linking model name (character string) to a list describing the prior distribution of each parameter to be estimated
#' @param ss_obs the observed summary statistics
#' @param max_number_of_gen the maximum number of generations to be performed
#' @param nb_acc_prtcl_per_gen the number of particles (per model, the total number corresponding to this number multiplied by the number of models) to be accepted during a generation before moving on to the next one
#' @param var_prtrbtn_krnl the standard deviation to be used for the perturbation kernel, if NA the empirical standard deviation will be used instead
#' @param var_prtrbtn_krnl_min the minimum standard deviation to be used for the perturbation kernel (must be greater than zero)
#' @param scaling_sd_prtrbtn_krnl scaling parameter for increasing or decreasing the value of the standard deviation used for the perturbation kernel
#' @param model_jump_prob probability of changing model when creating a new particle
#' @param nb_threshold number of thresholds used, corresponding to the number of distances returned by the function(s) stipulated in the model list (model_list)
#' @param new_threshold_quantile the quantile to be used to decrease the threshold(s) during iterations
#' @param distance_threshold_min the minimum threshold below which the procedure stops (in order to prevent excessively long computations)
#' @param max_attempts the maximum number of particles to be tested during an iteration, beyond which the procedure stops (in order to prevent excessively long computations)
#' @param acceptance_rate_min the acceptance rate below which the procedure stops (in order to prevent excessively long computations)
#' @param use_lhs_for_first_iter whether or not to use a LHS sampling method for the first iteration
#' @param experiment_folderpath the folder in which to carry out the estimation procedure and save the results
#' @param on_cluster whether or not the procedure is run on a computation cluster
#' @param cluster_type cluster type used (sge and slurm currently supported)
#' @param slurm_script_template script used to launch jobs on a slurm cluster
#' @param sge_script_template script used to launch jobs on a sge cluster
#' @param max_concurrent_jobs maximum number of jobs/tasks run in parallel
#' @param abc_user_param_file_path an R file containing the algorithm's parameters (usage not recommended, included in this version for reasons of compatibility with the procedure in script form used in some projects)
#' @param previous_gens an object (dataframe) containing previous results (set of iterations), in order to start from the last iteration performed
#' @param previous_epsilons an object (dataframe) containing previous results (set of thresholds), in order to start from the last iteration performed
#' @param verbose whether or not to display specific information
#'
#' @return a list containing two dataframes corresponding to (1) the particles accepted and (2) the thresholds used, during the successive iterations
#' @export
#' @include createLHSfromPrior.R defineNextThreshold.R setEmpiricalSD.R subjob.R saveEnvir.R
#'
#' @examples
#' # library(BRREWABC)
#' #
#' # # model definition
#' # compute_dist = function(x, ss_obs){
#' #     ss_sim = c( x[["alpha"]] + x[["beta"]] + rnorm(1,0,0.1),
#' #            x[["alpha"]] * x[["beta"]] + rnorm(1,0,0.1) ) # a very simple toy model
#' #     dist = sum((ss_sim-ss_obs)^2)
#' #     return(c(dist))
#' # }
#' #
#' # MODEL_LIST <- list("m1" = compute_dist)
#' # PRIOR_DIST <- list("m1" = list(c('alpha', 'unif', 0, 4), c('beta', 'unif', 0, 1)))
#' #
#' # # create a reference trajectory
#' # sum_stat_obs = c(2.0,0.75)
#' #
#' # # run abc smc procedure
#' # res = abcsmc( model_def = "test_brrewabc_models.R", model_list = MODEL_LIST, prior_dist = PRIOR_DIST, ss_obs = sum_stat_obs, max_number_of_gen = 1, nb_acc_prtcl_per_gen = 10, new_threshold_quantile = 0.8, experiment_folderpath = "", max_concurrent_jobs = 1, verbose = TRUE )
#' #
#' # # get results and plots
#' # all_accepted_particles = res$particles
#' # all_thresholds = res$thresholds
#' # plot_abcsmc_res(data = all_accepted_particles, prior = PRIOR_DIST, filename = "pairplot_all.png", colorpal = "GnBu")
#' # plot_densityridges(data = all_accepted_particles, prior = PRIOR_DIST, filename = "densityridges.png", colorpal = "GnBu")
#' # plot_thresholds(data = all_thresholds, nb_threshold = 1, filename = "thresholds.png", colorpal = "GnBu")
abcsmc <- function( model_def = NULL, # required
                    model_list = list(), # required
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
# THE FOLLOWING SECTION SHOULD NOT BE MODIFIED
#SBATCH --job-name=job-array_%%A_%%a   # nom du job
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --hint=nomultithread
#SBATCH --time=24:00:00
#SBATCH --array=%s-%s%%%d
output_fpath=$s
error_fpath=$s
#SBATCH --output=$output_fpath/output_%%A_%%a.out
#SBATCH --error=$error_fpath/error_%%A_%%a.out
mkdir -p $output_fpath
mkdir -p $error_fpath
Rscript $s $SLURM_ARRAY_TASK_ID
', # TODO : queue selection via a function argument
                    sge_script_template = '#!/bin/bash
#$ -S /bin/bash
#$ -N subjob_abcsmc_prlll
#$ -q "short.q|long.q"
# THE FOLLOWING SECTION SHOULD NOT BE MODIFIED
#$ -cwd
#$ -V
#$ -t %s-%s
#$ -tc %d
#$ -o /dev/null
#$ -e /dev/null
output_fpath=$s
error_fpath=$s
mkdir -p $output_fpath
mkdir -p $error_fpath
Rscript $s $SGE_TASK_ID >$output_fpath/subjob.${SGE_TASK_ID}.out 2>$error_fpath/subjob.${SGE_TASK_ID}.err
', # TODO : queue selection via a function argument
                    max_concurrent_jobs = 1,
                    abc_user_param_file_path = NULL,
                    previous_gens = NA,
                    previous_epsilons = NA,
                    verbose = FALSE ) {

    # first, load user param file if exist
    if(!is.null(abc_user_param_file_path)) {
        source(abc_user_param_file_path)
    }

    #

    tmp_folder_path <- "tmp"
    results_folder_path <- "res"
    if (experiment_folderpath != "") {
        tmp_folder_path <- file.path(experiment_folderpath, tmp_folder_path)
        results_folder_path <- file.path(experiment_folderpath, results_folder_path)
    }

    tmp_accepted_particles_filepath <- file.path(tmp_folder_path, "tmp_current_gen_accepted_particles")
    tmp_all_tested_particles_filepath <- file.path(tmp_folder_path, "tmp_current_gen_all_tested_particles")
    tmp_local_task_std_out <- file.path(tmp_folder_path, "std_out")
    tmp_local_task_std_err <- file.path(tmp_folder_path, "std_err")
    tmp_current_abc_state <- file.path(tmp_folder_path, "currentABCState.RData")

    if (verbose) {
        print(tmp_current_abc_state)
    }

    results_folder_path_CSV <- file.path(results_folder_path, 'csv')
    results_folder_path_FIGS <- file.path(results_folder_path, 'figs')
    accepted_particles_filepath <- file.path(results_folder_path_CSV, "all_accepted_particles.csv")
    last_accepted_particles_filepath <- file.path(results_folder_path_CSV, "last_accepted_particles.csv")
    all_particles_filepath <- file.path(results_folder_path_CSV, "all_particles.csv")
    thresholds_filepath <- file.path(results_folder_path_CSV, "thresholds.csv")
    lhs_used_for_first_gen_filepath <- file.path(results_folder_path_CSV, "lhs_for_first_gen.csv")

    subjob_script_name <- "brrewabc_subjob.R"
    subjob_script_path <- file.path(tmp_folder_path, subjob_script_name)

    #

    # TODO : use a main folder with default name based on the datehoursminutesseconds if not specified
    for (folder_path in c(tmp_folder_path, results_folder_path, results_folder_path_CSV, results_folder_path_FIGS)) {
        if (verbose) {cat(paste0("Check folder_path for : ", folder_path, "\n"))}
        # # Extract the directory part of the file path
        # folder_path <- dirname(file_path)
        # # Check if the last character is "/"
        # path_last_character <- substr(file_path, nchar(file_path), nchar(file_path))
        # if (path_last_character == "/") {
        #   folder_path <- file_path
        # }
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

    subjob_script_content <- "library(BRREWABC\nargs <- commandArgs(trailingOnly = TRUE)\nid <- as.integer(args[1])\nsubjob(job_id = id, path_to_ABC_state = '%s')\n"
    subjob_script <- sprintf(subjob_script_content, tmp_current_abc_state)
    writeLines(subjob_script, subjob_script_path)


    dist_names <- paste0("dist", as.character(seq(1,nb_threshold,1)))
    model_names <- names(model_list)
    param_names <- unique(Reduce(c, sapply(prior_dist, function(x) sapply(x, `[[`, 1))))
    column_names <- c("gen", "model", param_names, "pWeight", dist_names)

    # define the total number of particles to accept before next gen, that will be
    # used as a upper limit for the number of simulation to run (avoid a while loop
    # without stopping criterion)
    nb_acc_prtcl_before_next_gen <- nb_acc_prtcl_per_gen*length(model_names)

    # define the vector of thresholds
    epsilon <- c()

    # define the vector of sd to use in the perturbation kernel
    empirical_sd = list()
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
    #
    #
    #


    #
    gen = 1
    epsilon_not_improved = 0

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
        all_acc_particles = previous_gens
        epsilons = previous_epsilons
        gen = max(all_acc_particles$gen)
        # previous_acc_particles : set based on last in all_acc_particles
        previous_acc_particles = as.data.frame(all_acc_particles[all_acc_particles$gen == gen,])
        # current_acc_particles  # nothing to do
        # epsilon
        epsilon = defineNextThreshold(previous_acc_particles, nb_threshold, dist_names, new_threshold_quantile)
        # calculate the empirical sd to be used in the perturbation kernel
        empirical_sd = setEmpiricalSd(previous_acc_particles, model_names, prior_dist, var_prtrbtn_krnl, scaling_sd_prtrbtn_krnl, var_prtrbtn_krnl_min)
        # set the generation number
        gen = gen + 1
    }
    #
    while (gen <= max_number_of_gen) {
        # rep = 1
        nb_accepted = 0
        totattempts = 0
        if (verbose) {
            cat("gen", gen, '\n')
            cat("threshold:", epsilon, '\n')
            cat("prtrbtn_krnl_sd:", unlist(empirical_sd), '\n')
        }
        #
        tmp_acc_prtcls = stats::setNames(data.frame(matrix(ncol = length(column_names), nrow = 0), stringsAsFactors=FALSE), column_names)
        tmp_all_prtcls = stats::setNames(data.frame(matrix(ncol = length(column_names), nrow = 0), stringsAsFactors=FALSE), column_names)
        utils::write.csv(tmp_acc_prtcls, paste0(tmp_accepted_particles_filepath,"_",gen,".csv"), row.names=FALSE, quote=FALSE)
        utils::write.csv(tmp_all_prtcls, paste0(tmp_all_tested_particles_filepath,"_",gen,".csv"), row.names=FALSE, quote=FALSE)

        # print(ls()) # DEBUG
        var_to_save = c( "model_def", "model_list", "prior_dist", "ss_obs", "model_jump_prob", "use_lhs_for_first_iter", "max_concurrent_jobs", "tmp_accepted_particles_filepath", "tmp_all_tested_particles_filepath", "dist_names", "model_names", "param_names", "column_names", "nb_acc_prtcl_before_next_gen", "epsilon", "empirical_sd", "lhs_first_gen", "gen", "previous_acc_particles")
        # saveEnvir(var_to_save, tmp_current_abc_state) # TODO : using a function is not working, need to fix this
        do.call("save", c(var_to_save, list(file = tmp_current_abc_state)))
        # save.image(tmp_current_abc_state) # DEBUG

        #
        cluster_job_id = NA
        local_processes = c()
        if (on_cluster) {
            # submit a job array with as many job as nb_acc_prtcl_before_next_gen (for gen == 1, need at one job per each particle)
            cluster_script = ""
            # Replace placeholders in the template with actual parameter values
            nbjobs = max_concurrent_jobs
            if (gen == 1) {nbjobs = min(nb_acc_prtcl_before_next_gen, max_concurrent_jobs)}
            if (cluster_type == "slurm") {
                cluster_script <- sprintf(slurm_script_template, 1, nbjobs, max_concurrent_jobs, tmp_local_task_std_out, tmp_local_task_std_err, subjob_script_path)
            } else if (cluster_type == "sge") {
                cluster_script <- sprintf(sge_script_template, 1, nbjobs, max_concurrent_jobs, tmp_local_task_std_out, tmp_local_task_std_err, subjob_script_path)
            } else {
                stop("Cluster type not supported in the current version")
            }
            # Write the modified cluster job script to a file
            writeLines(cluster_script, "abc_smc_array_job.sh")
            # Submit the array job and capture the job ID
            if (cluster_type == "slurm") {
                cluster_job_id <- system('sbatch abc_smc_array_job.sh | awk \'{print $4}\'', intern = TRUE)
                print(paste("Submitted Slurm array job with ID:", cluster_job_id))
            } else if (cluster_type == "sge") {
                cluster_job_id <- system('qsub abc_smc_array_job.sh | awk -F "." \'{print $1}\' | awk \'{print $3}\'', intern = TRUE)
                print(paste("Submitted SGE array job with ID:", cluster_job_id))
            } else {
                stop("Cluster type not supported in the current version")
            }
        } else {
            # Parallel task on local machine
            for (task_index in 1:max_concurrent_jobs) {
                # Launch external script
                local_process = callr::r_bg(
                # local_process = callr::r(
                    # func = subjob,
                    func = function() {
                        library(BRREWABC)  # Load your package
                        subjob(job_id = task_index, path_to_ABC_state = tmp_current_abc_state)  # Call your function
                    },
                    # args = list(job_id = task_index, path_to_ABC_state = tmp_current_abc_state),
                    stdout = paste0(tmp_local_task_std_out,"abc_smc_task_",task_index,".out"),
                    stderr = paste0(tmp_local_task_std_err,"abc_smc_task_",task_index,".err"),
                    package = TRUE
                )
                local_processes <- c(local_processes, local_process)
            }
        }
        #
        # Create a progress bar
        pb <- progress::progress_bar$new(format = "gen :gen [:bar] :percent (ar: :accrate | :nbattempt) | eta: :eta (:elapsed)", clear = FALSE, total = nb_acc_prtcl_before_next_gen)
        #
        current_iter_broke = FALSE
        while (nb_accepted < nb_acc_prtcl_before_next_gen) { # repeat until N particles accepted
            tmp_accepted_particles = utils::read.csv(paste0(tmp_accepted_particles_filepath,"_",gen,".csv"))
            tmp_all_tested_particles = utils::read.csv(paste0(tmp_all_tested_particles_filepath,"_",gen,".csv"))
            nb_accepted = nrow(tmp_accepted_particles)
            totattempts = nrow(tmp_all_tested_particles)
            current_acc_rate = nb_accepted/totattempts
            # Increment the progress bar
            pb$update(min(nb_accepted/nb_acc_prtcl_before_next_gen, 1.00), tokens = list(gen = gen, accrate = format(round(current_acc_rate,digits=3),nsmall=3), nbattempt = totattempts))
            #
            if (totattempts > max_attempts) {
                message('\n', "The maximum number of attempts to accept the N particles has been reached!", '\n')
                current_iter_broke = TRUE
                break
            }
            if ((totattempts >= (nb_acc_prtcl_per_gen*length(model_names))) & (current_acc_rate < acceptance_rate_min)) {
                message('\n', "The acceptance rate has become too low (< specified acceptance_rate_min), the algorithm stops!")
                current_iter_broke = TRUE
                break
            }
        }
        # cancel subjob on cluster
        if (gen == 1) {
            # nothing to do
        } else {
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
        }
        # Close the progress bar
        # pb$terminate()
        # cat('\n')
        #
        if (current_iter_broke == TRUE) {break}
        #
        current_acc_particles = utils::read.csv(paste0(tmp_accepted_particles_filepath,"_",gen,".csv"))
        current_acc_particles = current_acc_particles[1:min(nrow(current_acc_particles),nb_acc_prtcl_before_next_gen),] # keep only the number of particle needed # TODO : imporve comment
        # normalise the weights
        for (mm in model_names) {
            if (mm %in% unique(current_acc_particles$model)) {
        	   current_acc_particles[current_acc_particles$model == mm,"pWeight"]= current_acc_particles[current_acc_particles$model == mm,"pWeight"]/ sum(current_acc_particles[current_acc_particles$model == mm,"pWeight"])
            }
        }
        utils::write.csv(current_acc_particles, last_accepted_particles_filepath, row.names=FALSE, quote=FALSE)
        # switch current and previous particles
        previous_acc_particles = current_acc_particles
        # save current particles
        all_acc_particles = rbind(all_acc_particles, current_acc_particles)
        # clear current particles table for next generation
        current_acc_particles = current_acc_particles[c(),]
        #
        if (gen > 1) {
            # save current epsilon if generation completed
            epsilons[nrow(epsilons) + 1,] = c(gen, epsilon, current_acc_rate)
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
        if ((gen > 1) & (all(epsilon <= distance_threshold_min))) {
            message("The distance threshold(s) (epsilon(s)) fall(s) below the predetermined min value!")
            print(epsilon)
            break
        }
        # calculate the empirical sd to be used in the perturbation kernel
        empirical_sd = setEmpiricalSd(previous_acc_particles, model_names, prior_dist, var_prtrbtn_krnl, scaling_sd_prtrbtn_krnl, var_prtrbtn_krnl_min)
        # define the next threshold
        next_epsilon = defineNextThreshold(previous_acc_particles, nb_threshold, dist_names, new_threshold_quantile)
        if ((gen > 1) & (all(next_epsilon == epsilon))) {
            epsilon_not_improved = epsilon_not_improved + 1
            if (epsilon_not_improved >= 2) {
                message("The distance threshold(s) (epsilon(s)) has not been improved over the last two iterations!")
                print(epsilons)
                break
            }
        }
        epsilon = next_epsilon
        #
        gen = gen + 1
        #
        # if (verbose) {cat('\n-\n')}
        if (verbose) {cat('-\n')}
    }
    if (verbose) {
        cat("Experiment done!", '\n')
        # cat("all thresholds:", '\n')
        # print(epsilons)
        # print(summary(previous_acc_particles[,c("model", param_names)]))
    }
    # cleaning
    if (on_cluster) {
        unlink(subjob_script_path)
        unlink("abc_smc_array_job.sh")
    } else {
        # for (task_index in 1:max_concurrent_jobs) {
        #     unlink(paste0("abc_smc_task_", task_index, ".sh"))
        # }
    }
    # rm tmp ?
    unlink(tmp_folder_path, recursive = TRUE)
    #
    return(list("particles" = all_acc_particles, "thresholds" = epsilons))
}
