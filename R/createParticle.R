
# #' Create a particle
# #'
# #' @param gen
# #' @param rep
# #' @param lhs_for_first_iter
# #' @param previous_acc_particles
# #' @param empirical_sd
# #' @param use_lhs_for_first_iter
# #' @param prior_dist
# #' @param model_names
# #' @param model_jump_prob
# #'
# #' @return
# # ' @export
# #'
# #' @examples
createParticle <- function(gen,
                           rep,
                           lhs_for_first_iter,
                           previous_acc_particles,
                           empirical_sd,
                           use_lhs_for_first_iter,
                           prior_dist,
                           model_names,
                           model_jump_prob) {
  the_particle <- list()
  if (gen == 1) {
    # sample from prior of get from lhs
    if (use_lhs_for_first_iter) {
      the_particle <- as.list(lhs_for_first_iter[rep, ])
    } else {
      the_particle[["model"]] <- sample(model_names, 1)
      for (pp in prior_dist[[the_particle[["model"]]]]) {
        if (pp[2] == "unif") {
          val <- stats::runif(1, min = as.double(pp[3]), max = as.double(pp[4]))
          the_particle[[pp[1]]] <- val
        } else {
          stop(p[2] + "- type of distribution not supported in the current version")
        }
      }
    }
  } else {
    # sample a model from previous population
    the_particle[["model"]] <- sample(previous_acc_particles$model, 1)
    # perturb the model
    if (stats::rbinom(1, 1, model_jump_prob)) {
      the_particle[["model"]] <- sample(unique(previous_acc_particles$model), 1)
    }
    # select all particle corresponding to the chosen model, in the previous population
    previous_acc_particles_subpop <- previous_acc_particles[previous_acc_particles$model == the_particle[["model"]], ]
    # sample from previous population with associated weight
    samp_idx <- sample(seq_len(nrow(previous_acc_particles_subpop)), 1, prob = previous_acc_particles_subpop$pWeight)
    sampled_particle <- previous_acc_particles_subpop[samp_idx, ]
    the_particle <- sampled_particle[c("model", sapply(prior_dist[[the_particle[["model"]]]], "[[", 1))]
    # perturb the particle
    for (pp in prior_dist[[the_particle[["model"]]]]) {
      if (pp[2] == "unif") {
        while (1 != 0) {
          #print(empirical_sd)
          prtrbd_val <- stats::rnorm(1, mean = the_particle[[pp[1]]],
                                     sd = empirical_sd[[the_particle[["model"]]]][[pp[1]]])
          # check if new value contained in initial prior
          if (stats::dunif(prtrbd_val, min = as.double(pp[3]), max = as.double(pp[4])) > 0) {
            the_particle[[pp[1]]] <- prtrbd_val
            break
          }
        }
      } else {
        stop(p[2] + " - type of distribution not supported in the current version")
      }
    }
  }
  return(the_particle)
}
