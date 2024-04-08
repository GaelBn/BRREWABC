
# #' Define next threshold(s)
# #'
# #' @param acc_particles
# #' @param nb_threshold
# #' @param distnames
# #' @param new_threshold_quantile
# #'
# #' @return
# # ' @export
# #'
# #' @examples
defineNextThreshold <- function(acc_particles, nb_threshold, distnames, new_threshold_quantile){
    the_epsilon = c()
    if (nb_threshold > 1) {
        the_epsilon = unname(apply( acc_particles[,distnames], 2, stats::quantile, probs = new_threshold_quantile, na.rm = TRUE))
    } else {
        the_epsilon = c(unname(stats::quantile(acc_particles$dist, new_threshold_quantile)))
    }
    return(the_epsilon)
}
