
# #' Run a subtask of the ABC-SMC
# #'
# #' @param job_id
# #' @param path_to_ABC_state
# #'
# #' @return
# # ' @export
# #' @include createParticle.R computeWeight.R saveEnvir.R
# #'
# #' @examples
subjob <- function( job_id,
                    path_to_ABC_state ) {

  # library(BRREWABC)

  # loadEnvir(path_to_ABC_state) # TODO : using a  function is not working, need to fix this
  load(path_to_ABC_state)
  source(model_def)

  #           _                      _       _     _
  #  ___  ___| |_   __   ____ _ _ __(_) __ _| |__ | | ___  ___
  # / __|/ _ \ __|  \ \ / / _` | '__| |/ _` | '_ \| |/ _ \/ __|
  # \__ \  __/ |_    \ V / (_| | |  | | (_| | |_) | |  __/\__ \
  # |___/\___|\__|    \_/ \__,_|_|  |_|\__,_|_.__/|_|\___||___/

  # set some variables depending on the generation id
  current_gen_acc_prtcls = 1
  current_job_nb_acc_prtcl_before_next_gen = nb_acc_prtcl_before_next_gen
  if (gen == 1) {
    current_gen_acc_prtcls = 1 + (job_id-1)*(nb_acc_prtcl_before_next_gen%/%min(max_concurrent_jobs,nb_acc_prtcl_before_next_gen))
    if (job_id < min(max_concurrent_jobs,nb_acc_prtcl_before_next_gen)) {
      current_job_nb_acc_prtcl_before_next_gen = current_gen_acc_prtcls + (nb_acc_prtcl_before_next_gen%/%min(max_concurrent_jobs,nb_acc_prtcl_before_next_gen)) - 1
    }
  }

  print(paste(job_id, current_gen_acc_prtcls, current_job_nb_acc_prtcl_before_next_gen, sep=" - "))


  #  _                 _        _       _
  # | | ___   __ _  __| |    __| | __ _| |_ __ _
  # | |/ _ \ / _` |/ _` |   / _` |/ _` | __/ _` |
  # | | (_) | (_| | (_| |  | (_| | (_| | || (_| |
  # |_|\___/ \__,_|\__,_|   \__,_|\__,_|\__\__,_|


  #  _            _                       _   _      _
  # | |_ ___  ___| |_    _ __   __ _ _ __| |_(_) ___| | ___  ___
  # | __/ _ \/ __| __|  | '_ \ / _` | '__| __| |/ __| |/ _ \/ __|
  # | ||  __/\__ \ |_   | |_) | (_| | |  | |_| | (__| |  __/\__ \
  #  \__\___||___/\__|  | .__/ \__,_|_|   \__|_|\___|_|\___||___/
  #                     |_|

  while (current_gen_acc_prtcls <= current_job_nb_acc_prtcl_before_next_gen) {
    # create the particle
    proposed_particle = createParticle(gen, current_gen_acc_prtcls, lhs_first_gen, previous_acc_particles, empirical_sd, use_lhs_for_first_iter, prior_dist, model_names, model_jump_prob)
    proposed_particle[["job_id"]] = job_id # if job_id need to be used in the model function
    # simulate and compute the distance
    dist = model_list[[proposed_particle[["model"]]]](proposed_particle, ss_obs)

    # compute the weight
    pWeight = NA
    if (gen == 1) {
        # compute the weight
        pWeight = 1
    } else {
        # accept particle if d <= epsilon, and so compute the weight
        if (all(dist <= epsilon)) {
            # compute the weight
            pWeight = computeWeight(proposed_particle, previous_acc_particles, empirical_sd, prior_dist)
        }
    }

    # save particle in the shared table
    if ((gen == 1) || (all(dist <= epsilon))) {
      # build the new row
      new.row = c(list(gen=gen, pWeight=pWeight), proposed_particle, stats::setNames(as.list(dist), dist_names))
      new.row = data.frame(new.row)
      missing_columns = setdiff(column_names, colnames(new.row)) # add the missing columns to the new row with empty values (or NA)
      new.row[missing_columns] = NA  # You can also define other default values if required
      new.row = new.row[column_names] # sort the columns to keep them in the right order
      # put filelock::lock on resfile.csv.lck
      lck = filelock::lock(paste0(tmp_accepted_particles_filepath,"_",gen,".csv.lck"))
      # Write the new line to the CSV file without reading it first
      utils::write.table(new.row, file = paste0(tmp_accepted_particles_filepath,"_",gen,".csv"), sep = ",", append = TRUE, col.names = FALSE, row.names = FALSE)
      # remove filelock::lock on resfile.csv.lck
      filelock::unlock(lck)
      # increment the number of accepted particle
      current_gen_acc_prtcls = current_gen_acc_prtcls + 1
    }

    # in any case, add the tested particle in a file
    # build the new row
    new.row = c(list(gen=gen, pWeight=pWeight), proposed_particle, stats::setNames(as.list(dist), dist_names))
    new.row = data.frame(new.row)
    missing_columns = setdiff(column_names, colnames(new.row)) # add the missing columns to the new row with empty values (or NA)
    new.row[missing_columns] = NA  # You can also define other default values if required
    new.row = new.row[column_names] # sort the columns to keep them in the right order
    # put filelock::lock on resfile.csv.lck
    lck = filelock::lock(paste0(tmp_all_tested_particles_filepath,"_",gen,".csv.lck"))
    # Write the new line to the CSV file without reading it first
    utils::write.table(new.row, file = paste0(tmp_all_tested_particles_filepath,"_",gen,".csv"), sep = ",", append = TRUE, col.names = FALSE, row.names = FALSE)
    # remove filelock::lock on resfile.csv.lck
    filelock::unlock(lck)

  }
}
