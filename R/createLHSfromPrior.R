
# #' Create LHS from prior distribution(s)
# #'
# #' @param model_names
# #' @param prior_dist
# #' @param nb_acc_prtcl_per_gen
# #'
# #' @return
# # ' @export
# #'
# #' @examples
createLHSfromPrior <- function(model_names, prior_dist, nb_acc_prtcl_per_gen) {
  all_lhs_for_first_iter <- list()
  for (mm in model_names) {
    tmp_lhs <- as.data.frame(lhs::randomLHS(nb_acc_prtcl_per_gen, length(prior_dist[[mm]])))
    the_param_names <- sapply(prior_dist[[mm]], `[[`, 1)
    colnames(tmp_lhs) <- the_param_names
    tmp_lhs["model"] <- mm
    for (pp in prior_dist[[mm]]) {
      tmp_lhs[, pp[1]] = scales::rescale(tmp_lhs[, pp[1]], to = c(as.double(pp[3]), as.double(pp[4])), from = c(0.0, 1.0))
    }
    all_lhs_for_first_iter[[mm]] <- tmp_lhs
  }
  the_lhs_for_first_iter <- dplyr::bind_rows(all_lhs_for_first_iter)
  return(the_lhs_for_first_iter)
}
