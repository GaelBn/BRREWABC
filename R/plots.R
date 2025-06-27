
#' Plot abcsmc results : pairplot and model selection
#'
#' @param data a dataframe containing the estimation results (set of particles
#' accepted during iterations)
#' @param prior a list linking model name (character string) to a list
#' describing the prior distribution of each parameter estimated (the same as
#' used for the abcsmc function)
#' @param filename the file name to be used to save the plots (the extension
#' defines the format: pdf or png)
#' @param figtitle the figure title
#' @param colorpal a palette name as used in the RColorBrewer package
#' @param iter if not NA, the id (number) of the iteration to be plotted
#'
#' @return one or several plots in pdf or png format
#' @export
#'
#' @examples
#' # see the abcsmc function help for details on how to plot results
plot_abcsmc_res <- function(data,
                            prior,
                            filename = "pairplot.png",
                            figtitle = "ABC smc results",
                            colorpal = "Greys",
                            iter = NA) {
  # message("The plot may fail if, for one or more iterations,
  # only one particle has been selected for one of the models.")
  tmp <- data
  colorshift <- 2 # because the first colours is often too light
  nb_gen = max(tmp$gen)
  nb_cols <- nb_gen + colorshift
  mycolors <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(8, colorpal))(nb_cols)
  mycolors <- mycolors[seq(-1, -colorshift, -1)]
  gen_threshold <- 15
  # TODO : add warning if nb gen > gen_threshold
  if((nb_gen > gen_threshold) || (length(iter) > gen_threshold)) {
    print("Number of generations exceed the threshold (15) allowed by ggpairs, it may cause long processing times. You may (re)define the iter argument to choose which generations to plot.")
  }
  # if((nb_gen > gen_threshold) || (length(iter) > gen_threshold)) {
  #   print("Number of generations exceed the threshold (15) allowed by ggpairs, a reduced number of generations have been plotted. You may (re)define the iter argument to choose which generations to plot.")
  #   modulo_nb_gen <- ceiling(nb_gen / gen_threshold)
  #   iter <- seq(nb_gen, 1, -modulo_nb_gen)
  # }
  if (all(!is.na(iter))) {
    tmp <- data[data$gen %in% iter, ]
    mycolors <- mycolors[iter]
  }

  tmp$gen <- as.factor(tmp$gen)
  # Get the last three characters
  last_three <- substr(filename, nchar(filename) - 2, nchar(filename))
  if (last_three == "png") {
    # Save plot as PNG
    print("Plot saved as 'png'.")
  } else if (last_three == "pdf") {
    # Save plot as PDF
    grDevices::pdf(filename, width = 9, height = 9)
    print("Plot saved as 'pdf'.")
  } else {
    stop("The specified file format is neither 'png' nor 'pdf'.")
  }
  # param_names = unique(Reduce(c, sapply(prior, function(x) sapply(x, `[[`, 1))))
  # pairplot = ggpairs(tmp[c(param_names, "gen")], aes(fill = gen, colour = gen, alpha=0.8),
  #     columns = param_names, lower = list(continuous = wrap("points", alpha = 0.5, size=1.5)),
  #     title = figtitle
  #  ) + scale_fill_manual(values = mycolors) + scale_colour_manual(values = mycolors) + theme_dark()
  # print(pairplot)
  if (length(unique(tmp$model)) >= 1) {
    if (length(unique(tmp$model)) > 1) {
      tmpm <- tmp[, c("gen", "model")]
      modelcolors <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(8, colorpal))(length(unique(tmpm$model)) + 1)
      modelcolors <- modelcolors[-1]
      stackedbarplot <- ggplot2::ggplot(data, ggplot2::aes(x = factor(gen), fill = model)) +
        ggplot2::geom_bar(position = "fill") +
        ggplot2::labs(x = "gen", y = "Count", fill = "model") +
        ggplot2::ggtitle(paste(figtitle, "model selection", sep = " - ")) +
        ggplot2::scale_fill_manual(values = modelcolors) +
        ggplot2::theme_minimal() +
        ggplot2::theme(legend.position = "bottom", panel.grid.minor = ggplot2::element_blank(), panel.grid.major = ggplot2::element_blank(), panel.grid.major.y = ggplot2::element_line(linewidth = .1, color = "black"))
      if (last_three == "png") {
        newfilename <- paste0(substr(filename, 1, nchar(filename) - 4),
                              "_modelprop", ".png")
        grDevices::png(newfilename, width = 9, height = 9, units = "in", res = 150)
      }
      print(stackedbarplot)
      if (last_three == "png") {
        grDevices::dev.off()
      }
    }
    for (mm in sort(unique(tmp$model))) {
      param_names <- sapply(prior[[mm]], "[[", 1)
      tmpm <- tmp[tmp$model == mm, ]
      current_fig_title = paste(figtitle, mm, sep = " - ")
      if (length(unique(tmp$model)) == 1) {
        current_fig_title = figtitle
      }
      pairplot <- GGally::ggpairs(tmpm[c(param_names, "gen")], upper = "blank", ggplot2::aes(fill = gen, colour = gen, alpha = 0.8), columns = param_names, lower = list(continuous = GGally::wrap("points", alpha = 0.5, size = 1.5)), title = current_fig_title, cardinality_threshold = NULL) +
        ggplot2::scale_fill_manual(values = mycolors) + ggplot2::scale_colour_manual(values = mycolors) + ggplot2::theme_minimal()
      if (last_three == "png") {
        newfilename <- paste0(substr(filename, 1, nchar(filename) - 4),
                              "_", mm, ".png")
        if (length(unique(tmp$model)) == 1) {
          newfilename <- filename
        }
        grDevices::png(newfilename, width = 9, height = 9,
                       units = "in", res = 150)
      }
      print(pairplot)
      if (last_three == "png") {
        grDevices::dev.off()
      }
    }
  }
  if (last_three == "pdf") {
    grDevices::dev.off()
  }
}


#' Plot abcsmc results : densityridges for estimated parameters
#'
#' @param data a dataframe containing the estimation results (set of particles
#' accepted during iterations)
#' @param prior a list linking model name (character string) to a list
#' describing the prior distribution of each parameter estimated (the same as
#' used for the abcsmc function)
#' @param filename the file name to be used to save the plots (the extension
#' defines the format: pdf or png)
#' @param figtitle the figure title
#' @param colorpal a palette name as used in the RColorBrewer package
#'
#' @return one or several plots in pdf or png format
#' @export
#'
#' @examples
#' # see the abcsmc function help for details on how to plot results
plot_densityridges <- function(data,
                               prior,
                               filename = "densityridges.png",
                               figtitle = "",
                               colorpal = "Greys") {
  tmp <- data
  colorshift <- 2  # because the first colours is often too light
  nb_cols <- max(tmp$gen) + 1 + colorshift
  mycolors <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(8, colorpal))(nb_cols)
  mycolors <- mycolors[seq(-1, -colorshift, -1)]
  tmp$gen <- as.factor(tmp$gen)
  # Get the last three characters
  last_three <- substr(filename, nchar(filename) - 2, nchar(filename))
  if (last_three == "png") {
    # Save plot as PNG
    print("Plot saved as 'png'.")
  } else if (last_three == "pdf") {
    # Save plot as PDF
    grDevices::pdf(filename, width = 6, height = 8)
    print("Plot saved as 'pdf'.")
  } else {
    stop("The specified file format is neither 'png' nor 'pdf'.")
  }
  for (mm in sort(unique(tmp$model))) {
    for (p in prior[[mm]]) {
      subtmp <- tmp[tmp$model == mm, c("gen", p[1])]
      names(subtmp)[names(subtmp) == p[1]] <- "param"
      sub_plot <- ggplot2::ggplot(subtmp, ggplot2::aes(x = param, y = gen, fill = gen, color = gen))+
        ggridges::geom_density_ridges(ggplot2::aes(), scale = 2.0, alpha = .6, bandwidth = (as.double(p[4]) - as.double(p[3])) / 100) +
        ggridges::theme_ridges() +
        ggplot2::scale_y_discrete(expand = ggplot2::expansion(add = c(0.3, 2.5))) +
        ggplot2::scale_x_continuous(limits = c(as.double(p[3]), as.double(p[4])), expand = c(0.01, 0)) +
        ggplot2::scale_fill_manual(values = mycolors) +
        ggplot2::scale_colour_manual(values = mycolors) +
        ggplot2::labs(x = p[1], y = "iteration") +
        # ggplot2::ggtitle(paste("Successive posterior distributions for parameter",p[1], sep=' '))
        ggplot2::theme_minimal() + ggplot2::theme(legend.position = "none")
      if (last_three == "png") {
          newfilename <- paste0(substr(filename, 1, nchar(filename) - 4), "_", mm, "_", p[1], ".png")
          if (length(unique(tmp$model)) == 1) {
            newfilename <- paste0(substr(filename, 1, nchar(filename) - 4), "_", p[1], ".png")
          }
          grDevices::png(newfilename, width = 6, height = 8, units = "in", res = 150)
      }
      print(sub_plot)
      if (last_three == "png") {
        grDevices::dev.off()
      }
    }
  }
  if (last_three == "pdf") {
    grDevices::dev.off()
  }
}


#' Plot abcsmc results : thresholds over iterations
#'
#' @param data a dataframe containing the estimation results (thresholds over
#' iterations)
#' @param nb_threshold number of thresholds used (see the abcsmc function help
#' for more details)
#' @param filename the file name to be used to save the plots (the extension
#' defines the format: pdf or png)
#' @param figtitle the figure title
#' @param colorpal a palette name as used in the RColorBrewer package
#'
#' @return one or several plots in pdf or png format
#' @export
#'
#' @examples
#' # see the abcsmc function help for details on how to plot results
plot_thresholds <- function(data,
                            nb_threshold = 1,
                            filename = "thresholds.png",
                            figtitle = "",
                            colorpal = "Greys") {
  tmp <- data
  colorshift <- 2  # because the first colours is often too light
  nb_cols <- max(tmp$gen) + 1 + colorshift
  mycolors <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(8, colorpal))(nb_cols)
  mycolors <- mycolors[seq(-1, -colorshift, -1)]
  tmp$gen <- as.factor(tmp$gen)
  # Get the last three characters
  last_three <- substr(filename, nchar(filename) - 2, nchar(filename))
  if (last_three == "png") {
    # Save plot as PNG
    grDevices::png(filename, width = 8, height = 5, units = "in", res = 150)
    print("Plot saved as 'png'.")
  } else if (last_three == "pdf") {
    # Save plot as PDF
    grDevices::pdf(filename, width = 8, height = 5)
    print("Plot saved as 'pdf'.")
  } else {
    stop("The specified file format is neither 'png' nor 'pdf'.")
  }
  for (dd in paste0("dist", as.character(seq(1, nb_threshold, 1)))) {
    subtmp <- tmp[, c("gen", dd)]
    names(subtmp)[names(subtmp) == dd] <- "dist"
    sub_plot = ggplot2::ggplot(subtmp, ggplot2::aes(x = gen, y = dist, fill = gen, color = gen))+
      ggplot2::geom_bar(stat = "identity") +
      ggplot2::scale_fill_manual(values = mycolors) +
      ggplot2::scale_colour_manual(values = mycolors) +
      ggplot2::labs(y = "threshold", x = "iteration ID") +
      # ggplot2::ggtitle(paste("Successive posterior distributions for parameter",p[1], sep=' '))
      ggplot2::theme_minimal() + ggplot2::theme(legend.position = "none")
    if (last_three == "png") {
      newfilename <- paste0(substr(filename, 1, nchar(filename) - 4), "_", dd, ".png")
      if (nb_threshold == 1) {
        newfilename <- filename
      }
      grDevices::png(newfilename, width = 8, height = 5, units = "in", res = 150)
    }
    print(sub_plot)
    if (last_three == "png") {
      grDevices::dev.off()
    }
  }
  if (last_three == "pdf") {
    grDevices::dev.off()
  }
}


#' Plot abcsmc results : ESS (Effective Sample Size) over iterations
#'
#' @param data a dataframe containing the estimation results (thresholds over
#' iterations)
#' @param filename the file name to be used to save the plots (the extension
#' defines the format: pdf or png)
#' @param figtitle the figure title
#' @param colorpal a palette name as used in the RColorBrewer package
#'
#' @return a dataframe and one or several plots in pdf or png format
#' @export
#'
#' @examples
#' # see the abcsmc function help for details on how to plot results
plot_ess <- function(data,
                     filename = "ess.png",
                     figtitle = "",
                     colorpal = "Greys") {
  tmp <- data
  colorshift <- 2  # because the first colours is often too light
  nb_cols <- max(tmp$gen) + 1 + colorshift
  mycolors <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(8, colorpal))(nb_cols)
  mycolors <- mycolors[seq(-1, -colorshift, -1)]
  tmp$gen <- as.factor(tmp$gen)
  # Get the last three characters
  last_three <- substr(filename, nchar(filename) - 2, nchar(filename))
  if (last_three == "png") {
    # Save plot as PNG
    grDevices::png(filename, width = 8, height = 5, units = "in", res = 150)
    print("Plot saved as 'png'.")
  } else if (last_three == "pdf") {
    # Save plot as PDF
    grDevices::pdf(filename, width = 8, height = 5)
    print("Plot saved as 'pdf'.")
  } else {
    stop("The specified file format is neither 'png' nor 'pdf'.")
  }

  ess_per_gen <- aggregate(pWeight ~ gen, data = tmp, FUN = function(x) 1 / sum(x^2))
  colnames(ess_per_gen) <- c("gen", "ess")

  plot = ggplot2::ggplot(ess_per_gen, ggplot2::aes(x = gen, y = ess, fill = gen, color = gen))+
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::scale_fill_manual(values = mycolors) +
    ggplot2::scale_colour_manual(values = mycolors) +
    ggplot2::labs(y = "Effective Sample Size", x = "iteration ID") +
    ggplot2::theme_minimal() + ggplot2::theme(legend.position = "none")
  if (last_three == "png") {
    newfilename <- filename
    grDevices::png(newfilename, width = 8, height = 5, units = "in", res = 150)
  }
  print(plot)
  if (last_three == "png") {
    grDevices::dev.off()
  }

  if (last_three == "pdf") {
    grDevices::dev.off()
  }
  return(ess_per_gen)
}


#' Plot distribution of distances by model
#'
#' @param data a dataframe containing the estimation results (set of particles
#' accepted during iterations, including : gen, model, pWeight and one or more dist columns)
#' @param generation generation to plot; defaults to the maximum value in gen
#' @param filename the file name to save the plot (extension png or pdf determines format)
#' @param figtitle the figure title
#' @param colorpal a palette name as used in the RColorBrewer package
#' @param bins number of histogram bins
#' @param alpha histogram transparency level
#' @param adjust density curve smoothing adjustment
#'
#' @return a dataframe in long format of distances; the plot is saved to file
#' @export
#'
#' @examples
#' # plot_distances(data)
plot_distances <- function(data,
                           generation = NULL,
                           filename = "distances.png",
                           figtitle = "",
                           colorpal = "Set2",
                           bins = 20,
                           alpha = 0.6,
                           adjust = 1.5) {
  df <- data

  # select generation
  if (is.null(generation)) {
    generation <- max(df$gen, na.rm = TRUE)
  }
  df <- df[df$gen == generation, ]

  # identify distance columns and compute total distance
  dist_cols <- grep("^dist", names(df), value = TRUE)
  df$distTot <- rowSums(df[, dist_cols, drop = FALSE], na.rm = TRUE)
  dist_cols <- c(dist_cols, "distTot")

  # reshape to long format
  df_long <- tidyr::pivot_longer(df, cols = dist_cols, names_to = "distance", values_to = "value")

  # generate extended color palette based on number of models
  n_models <- length(unique(df_long$model))
  # ensure minimum palette size for RColorBrewer
  base_size <- max(3, min(8, n_models))
  base_cols <- RColorBrewer::brewer.pal(base_size, colorpal)
  # if fewer models than base_size, subset
  base_cols <- base_cols[seq_len(n_models)]
  mypal <- grDevices::colorRampPalette(base_cols)(n_models)

  # open graphics device
  ext <- tools::file_ext(filename)
  if (ext == "png") {
    grDevices::png(filename, width = 8, height = 5, units = "in", res = 150)
  } else if (ext == "pdf") {
    grDevices::pdf(filename, width = 8, height = 5)
  } else {
    stop("unsupported file format; use 'png' or 'pdf'")
  }

  # plot using updated aesthetics
  p <- ggplot2::ggplot(df_long, ggplot2::aes(x = value, fill = model)) +
    ggplot2::geom_histogram(ggplot2::aes(y = ggplot2::after_stat(density)), bins = bins, alpha = alpha, position = "identity", color = "black") +
    ggplot2::geom_density(alpha = alpha / 2, adjust = adjust) +
    ggplot2::facet_wrap(~ distance, scales = "free") +
    ggplot2::scale_fill_manual(values = mypal) +
    ggplot2::labs(x = "value", y = "density", title = figtitle) +
    ggplot2::theme_minimal()

  print(p)
  grDevices::dev.off()

  invisible(df_long)
}



#' Plot abcrejection results : pairplot
#'
#' @param data a dataframe containing the estimation results (set of particles
#' accepted during iterations)
#' @param prior a list linking model name (character string) to a list
#' describing the prior distribution of each parameter estimated (the same as
#' used for the abcsmc function)
#' @param thresholds if not NA, the threshold(s) (a value or a vector of
#' thresholds) to use to select acceptable particles. Should be use only if
#' abcrejection was run using thresholds = NA
#' @param filename the file name to be used to save the plots (the extension
#' defines the format: pdf or png)
#' @param figtitle the figure title
#' @param colorpal a palette name as used in the RColorBrewer package
#'
#' @return one or several plots in pdf or png format
#' @export
#'
#' @examples
#' # see the abcrejection function help for details on how to plot results
plot_abcrejection_res <- function(data,
                                  prior,
                                  thresholds = NA,
                                  filename = "pairplot.png",
                                  figtitle = "ABC rejection results",
                                  colorpal = "Greys") {
  # message("The plot may fail if, for one or more iterations,
  # only one particle has been selected for one of the models.")
  colorshift <- 2 # because the first colours is often too light
  nb_thresholds <- length(thresholds)
  nb_cols <- nb_thresholds + colorshift
  mycolors <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(8, colorpal))(nb_cols)
  mycolors <- mycolors[seq(-1, -colorshift, -1)]
  nb_thresholds_limit <- 15

  if(nb_thresholds > nb_thresholds_limit) {
    print("Number of thresholds exceed the limit (15) allowed by ggpairs, it may cause long processing times. You may (re)define the thresholds argument to choose which thresholds to plot.")
  }

  # Create a table based on the defined thresholds if not NA
  tmp = data.frame()
  if (all(!is.na(thresholds))) {
    list_of_dataframes <- list()
    for (tr in thresholds) {
      tmp_current_cutoff = data[data$dist1 <= tr,]
      tmp_current_cutoff$cutoff <- tr
      list_of_dataframes <- append(list_of_dataframes, list(tmp_current_cutoff))
    }
    tmp <- do.call(rbind, list_of_dataframes)
  } else {
    tmp <- data
    tmp$cutoff <- 0
  }
  tmp$cutoff <- factor(tmp$cutoff, levels = sort(unique(tmp$cutoff), decreasing = TRUE))

  # Get the last three characters
  last_three <- substr(filename, nchar(filename) - 2, nchar(filename))
  if (last_three == "png") {
    # Save plot as PNG
    print("Plot saved as '.png'.")
  } else if (last_three == "pdf") {
    # Save plot as PDF
    grDevices::pdf(filename, width = 9, height = 9)
    print("Plot saved as 'pdf'.")
  } else {
    stop("The specified file format is neither 'png' nor 'pdf'.")
  }
  # param_names = unique(Reduce(c, sapply(prior, function(x) sapply(x, `[[`, 1))))
  # pairplot = ggpairs(tmp[c(param_names, "cutoff")], aes(fill = cutoff, colour = cutoff, alpha=0.8),
  #     columns = param_names, lower = list(continuous = wrap("points", alpha = 0.5, size=1.5)),
  #     title = figtitle
  #  ) + scale_fill_manual(values = mycolors) + scale_colour_manual(values = mycolors) + theme_dark()
  # print(pairplot)
  if (length(unique(tmp$model)) >= 1) {
    if (length(unique(tmp$model)) > 1) {
      tmpm <- tmp[, c("cutoff", "model")]
      modelcolors <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(8, colorpal))(length(unique(tmpm$model)) + 1)
      modelcolors <- modelcolors[-1]
      stackedbarplot <- ggplot2::ggplot(data, ggplot2::aes(x = factor(cutoff), fill = model)) +
        ggplot2::geom_bar(position = "fill") +
        ggplot2::labs(x = "cutoff", y = "Count", fill = "model") +
        ggplot2::ggtitle(paste(figtitle, "model selection", sep = " - ")) +
        ggplot2::scale_fill_manual(values = modelcolors) +
        ggplot2::theme_minimal() +
        ggplot2::theme(legend.position = "bottom", panel.grid.minor = ggplot2::element_blank(), panel.grid.major = ggplot2::element_blank(), panel.grid.major.y = ggplot2::element_line(linewidth = .1, color = "black"))
      if (last_three == "png") {
        newfilename <- paste0(substr(filename, 1, nchar(filename) - 4),
                              "_modelprop", ".png")
        grDevices::png(newfilename, width = 9, height = 9, units = "in", res = 150)
      }
      print(stackedbarplot)
      if (last_three == "png") {
        grDevices::dev.off()
      }
    }
    for (mm in sort(unique(tmp$model))) {
      param_names <- sapply(prior[[mm]], "[[", 1)
      tmpm <- tmp[tmp$model == mm, ]
      current_fig_title = paste(figtitle, mm, sep = " - ")
      if (length(unique(tmp$model)) == 1) {
        current_fig_title = figtitle
      }
      pairplot <- GGally::ggpairs(tmpm[c(param_names, "cutoff")], upper = "blank", ggplot2::aes(fill = cutoff, colour = cutoff, alpha = 0.8), columns = param_names, lower = list(continuous = GGally::wrap("points", alpha = 0.5, size = 1.5)), title = current_fig_title, cardinality_threshold = NULL) +
        ggplot2::scale_fill_manual(values = mycolors) + ggplot2::scale_colour_manual(values = mycolors) + ggplot2::theme_minimal()
      if (last_three == "png") {
        newfilename <- paste0(substr(filename, 1, nchar(filename) - 4),
                              "_", mm, ".png")
        if (length(unique(tmp$model)) == 1) {
          newfilename <- filename
        }
        grDevices::png(newfilename, width = 9, height = 9,
                       units = "in", res = 150)
      }
      print(pairplot)
      if (last_three == "png") {
        grDevices::dev.off()
      }
    }
  }
  if (last_three == "pdf") {
    grDevices::dev.off()
  }
}
