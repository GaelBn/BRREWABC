
# #' Compute the empirical standard deviation to be used in perturbation kernels
# #'
# #' @param acc_particles
# #' @param model_names
# #' @param prior_dist
# #' @param var_prtrbtn_krnl
# #' @param scaling_sd_prtrbtn_krnl
# #' @param var_prtrbtn_krnl_min
# #'
# #' @return A vector of the empirical standard deviation to be used in perturbation kernels
# # ' @export
# #'
# #' @examples
setEmpiricalSd <- function(acc_particles,
                           model_names,
                           prior_dist,
                           var_prtrbtn_krnl,
                           scaling_sd_prtrbtn_krnl,
                           var_prtrbtn_krnl_min) {
  the_empirical_sd <- list()
  for (mm in model_names) {
    the_empirical_sd[[mm]] <- list()
    for (pp in prior_dist[[mm]]) {
      the_empirical_sd[[mm]][[pp[1]]] <- var_prtrbtn_krnl
    }
  }
  #
  if (is.na(var_prtrbtn_krnl)) {
    for (mm in model_names) {
      for (pp in prior_dist[[mm]]) {
        if (length(acc_particles[acc_particles$model == mm, pp[1]]) > 1) {
          the_empirical_sd[[mm]][[pp[1]]] <- stats::sd(acc_particles[acc_particles$model == mm, pp[1]]) * scaling_sd_prtrbtn_krnl[[pp[1]]]
        } else {
          the_empirical_sd[[mm]][[pp[1]]] <- var_prtrbtn_krnl_min
        }
      }
    }
  }
  return(the_empirical_sd)
}
