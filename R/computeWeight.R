
# #' Compute weight of a given particle
# #'
# #' @param proposed_particle
# #' @param previous_acc_particles
# #' @param empirical_sd
# #' @param prior_dist
# #'
# #' @return
# # ' @export
# #'
# #' @examples
computeWeight <- function(proposed_particle, previous_acc_particles, empirical_sd, prior_dist){
    model_used = proposed_particle[["model"]]
    pWeight = 0
    numerator = 1
    denominator = previous_acc_particles[previous_acc_particles$model == model_used,"pWeight"]
    for (pp in prior_dist[[model_used]]) {
        if (pp[2] == "unif") {
            numerator = numerator * stats::dunif(proposed_particle[[pp[1]]], min = as.double(pp[3]), max = as.double(pp[4]))
            denominator = denominator * stats::dnorm(proposed_particle[[pp[1]]], mean = previous_acc_particles[previous_acc_particles$model == model_used,pp[1]], sd = empirical_sd[[model_used]][[pp[1]]])
        } else {
            stop(pp[2] + " - type of distribution not supported in the current version")
        }
    }
    pWeight = numerator / sum(denominator)
    return(pWeight)
}
