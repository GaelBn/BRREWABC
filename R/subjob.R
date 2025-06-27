
#' Run a subtask of the ABC-SMC. Shouldn't have to be used by the user,
#' this function is visible so that it can be used on cluster by the main script
#'
#' @param job_id id of the current job
#' @param path_to_abc_state path to the .Rdata of the current experiment
#'
#' @return nothing, write results in a tmp file
#' @export
#' @include createParticle.R computeWeight.R saveEnvir.R
subjob_smc <- function(job_id,
                       path_to_abc_state) {

  # loadEnvir(path_to_abc_state) # TODO : not working, need to fix this
  load(path_to_abc_state)
  if (!is.null(model_def)) {
    source(model_def)
  }

  #           _                      _       _     _
  #  ___  ___| |_   __   ____ _ _ __(_) __ _| |__ | | ___  ___
  # / __|/ _ \ __|  \ \ / / _` | '__| |/ _` | '_ \| |/ _ \/ __|
  # \__ \  __/ |_    \ V / (_| | |  | | (_| | |_) | |  __/\__ \
  # |___/\___|\__|    \_/ \__,_|_|  |_|\__,_|_.__/|_|\___||___/

  # set some variables depending on the generation id
  nb_acc_prtcls <- 1
  current_job_nb_acc_prtcl_before_next_gen <- nb_acc_prtcl_before_next_gen
  if (gen == 1) {
    if (use_lhs_for_first_iter) {
      nbjobs = min(nb_acc_prtcl_before_next_gen, max_concurrent_jobs)
      q <- nb_acc_prtcl_before_next_gen %/% nbjobs    # base workload
      r <- nb_acc_prtcl_before_next_gen %%  nbjobs    # job with extra tasks
      if (job_id <= r) {
        start <- (job_id - 1) * (q + 1) + 1
        end   <- start + q          # q + 1 tasks
      } else {
        start <- r * (q + 1) + (job_id - r - 1) * q + 1
        end   <- start + q - 1      # q tasks
      }
      nb_acc_prtcls = start
      current_job_nb_acc_prtcl_before_next_gen = min(end, nb_acc_prtcl_before_next_gen)
    }
  }

  #  _            _                       _   _      _
  # | |_ ___  ___| |_    _ __   __ _ _ __| |_(_) ___| | ___  ___
  # | __/ _ \/ __| __|  | '_ \ / _` | '__| __| |/ __| |/ _ \/ __|
  # | ||  __/\__ \ |_   | |_) | (_| | |  | |_| | (__| |  __/\__ \
  #  \__\___||___/\__|  | .__/ \__,_|_|   \__|_|\___|_|\___||___/
  #                     |_|

  while (nb_acc_prtcls <= current_job_nb_acc_prtcl_before_next_gen) {
    # create the particle
    proposed_particle <- createParticleForSMC(gen, nb_acc_prtcls, lhs_first_gen, previous_acc_particles, empirical_sd, use_lhs_for_first_iter, prior_dist, model_names, model_jump_prob)
    proposed_particle[["job_id"]] <- job_id # if job_id need to be used in the model function
    # simulate and compute the distance
    dist <- model_list[[proposed_particle[["model"]]]](proposed_particle, ss_obs)

    # compute the weight
    pWeight <- NA
    if (gen == 1) {
      # compute the weight
      pWeight <- 1
    } else {
      # accept particle if d < epsilon, and so compute the weight
      if (all(dist < epsilon)) {
        # compute the weight
        pWeight <- computeWeight(proposed_particle, previous_acc_particles, empirical_sd, prior_dist)
      }
    }

    # save particle in the shared table
    if ((gen == 1) || (all(dist < epsilon))) {
      # build the new row
      new.row <- c(list(gen = gen, pWeight = pWeight), proposed_particle, stats::setNames(as.list(dist), dist_names))
      new.row <- data.frame(new.row)
      missing_columns <- setdiff(column_names, colnames(new.row)) # add the missing columns to the new row with empty values (or NA)
      new.row[missing_columns] <- NA  # You can also define other default values if required
      new.row <- new.row[column_names] # sort the columns to keep them in the right order
      # put filelock::lock on resfile.csv.lck
      lck <- filelock::lock(paste0(tmp_accepted_particles_filepath, "_", gen, ".csv.lck"))
      # Write the new line to the CSV file without reading it first
      utils::write.table(new.row, file = paste0(tmp_accepted_particles_filepath, "_", gen, ".csv"), sep = ",", append = TRUE, col.names = FALSE, row.names = FALSE)
      # remove filelock::lock on resfile.csv.lck
      filelock::unlock(lck)
      # increment the number of accepted particle
      nb_acc_prtcls <- nb_acc_prtcls + 1
    }

    # in any case, add the tested particle in a file
    # build the new row
    new.row <- c(list(gen = gen, pWeight = pWeight), proposed_particle, stats::setNames(as.list(dist), dist_names))
    new.row <- data.frame(new.row)
    missing_columns <- setdiff(column_names, colnames(new.row)) # add the missing columns to the new row with empty values (or NA)
    new.row[missing_columns] <- NA  # You can also define other default values if required
    new.row <- new.row[column_names] # sort the columns to keep them in the right order
    # put filelock::lock on resfile.csv.lck
    lck <- filelock::lock(paste0(tmp_all_tested_particles_filepath, "_", gen, ".csv.lck"))
    # Write the new line to the CSV file without reading it first
    utils::write.table(new.row, file = paste0(tmp_all_tested_particles_filepath, "_", gen, ".csv"), sep = ",", append = TRUE, col.names = FALSE, row.names = FALSE)
    # remove filelock::lock on resfile.csv.lck
    filelock::unlock(lck)

  }
}




#' Run a subtask of the ABC rejection Shouldn't have to be used by the user,
#' this function is visible so that it can be used on cluster by the main script
#'
#' @param job_id id of the current job
#' @param path_to_abc_state path to the .Rdata of the current experiment
#'
#' @return nothing, write results in a tmp file
#' @export
#' @include createParticle.R computeWeight.R saveEnvir.R
subjob_rejection <- function(job_id,
                   path_to_abc_state) {

  # loadEnvir(path_to_abc_state) # TODO : not working, need to fix this
  load(path_to_abc_state)
  if (!is.null(model_def)) {
    source(model_def)
  }

  #           _                      _       _     _
  #  ___  ___| |_   __   ____ _ _ __(_) __ _| |__ | | ___  ___
  # / __|/ _ \ __|  \ \ / / _` | '__| |/ _` | '_ \| |/ _ \/ __|
  # \__ \  __/ |_    \ V / (_| | |  | | (_| | |_) | |  __/\__ \
  # |___/\___|\__|    \_/ \__,_|_|  |_|\__,_|_.__/|_|\___||___/

  # set some variables depending on the generation id
  nb_acc_prtcls <- 0

  #  _            _                       _   _      _
  # | |_ ___  ___| |_    _ __   __ _ _ __| |_(_) ___| | ___  ___
  # | __/ _ \/ __| __|  | '_ \ / _` | '__| __| |/ __| |/ _ \/ __|
  # | ||  __/\__ \ |_   | |_) | (_| | |  | |_| | (__| |  __/\__ \
  #  \__\___||___/\__|  | .__/ \__,_|_|   \__|_|\___|_|\___||___/
  #                     |_|

  while (nb_acc_prtcls < tot_nb_acc_prtcl) {
    # create the particle
    proposed_particle <- createParticle(prior_dist, model_names)
    proposed_particle[["job_id"]] <- job_id # if job_id need to be used in the model function
    # simulate and compute the distance
    dist <- model_list[[proposed_particle[["model"]]]](proposed_particle, ss_obs)

    # save particle in the shared tables
    # build the new row
    new.row <- c(proposed_particle, stats::setNames(as.list(dist), dist_names))
    new.row <- data.frame(new.row)
    missing_columns <- setdiff(column_names, colnames(new.row)) # add the missing columns to the new row with empty values (or NA)
    new.row[missing_columns] <- NA  # You can also define other default values if required
    new.row <- new.row[column_names] # sort the columns to keep them in the right order

    # put filelock::lock on resfile.csv.lck
    lck <- filelock::lock(paste0(all_tested_particles_filepath, ".lck"))
    # Write the new line to the CSV file without reading it first
    utils::write.table(new.row, file = all_tested_particles_filepath, sep = ",", append = TRUE, col.names = FALSE, row.names = FALSE)
    # remove filelock::lock on resfile.csv.lck
    filelock::unlock(lck)

    if (all(!is.na(thresholds))) {
      if (all(dist <= thresholds)) {
        # put filelock::lock on resfile.csv.lck
        lck <- filelock::lock(paste0(accepted_particles_filepath, ".lck"))
        # Write the new line to the CSV file without reading it first
        utils::write.table(new.row, file = accepted_particles_filepath, sep = ",", append = TRUE, col.names = FALSE, row.names = FALSE)
        # remove filelock::lock on resfile.csv.lck
        filelock::unlock(lck)
        # increment the number of accepted particle
        nb_acc_prtcls <- nb_acc_prtcls + 1
      }
    }
  }
}
