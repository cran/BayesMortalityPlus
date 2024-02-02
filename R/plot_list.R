#' @name plot.list
#' @rdname plot.list
#'
#' @title Plot a set of life tables
#'
#' @description
#' Function that returns a log-scale 'ggplot' of the mortality graduation
#'  returned by hp(), dlm(), hp_close() or dlm_close() functions.
#'
#'
#' @param x List of objects of the following classes: `HP`, `DLM`, `ClosedHP` or `ClosedDLM`.
#' @param age Vector with the ages to plot the life tables.
#' @param Ex Vector with the exposures of the selected ages. Its length must be equal to the age vector. This argument is only necessary when plotting poisson and binomial HP models.
#' @param plotIC Logical. If 'TRUE'(default), plots the predictive intervals.
#' @param plotData Logical. If 'TRUE' (default), plots the data used in the modelling as dots.
#' @param labels Description of the curve (Optional).
#' @param colors Vector of colours of the curves (Optional).
#' @param linetype Vector with the line type of the curve. (Optional).
#' @param prob Coverage probability of the predictive intervals. Default is '0.95'.
#' @param ... Other arguments.
#'
#' @return A 'ggplot' object with fitted life tables.
#'
#' @examples
#' ## Importing mortality data from the USA available on the Human Mortality Database (HMD):
#' data(USA)
#'
#' ## Selecting the log mortality rate of the 1990 male population ranging from 0 to 100 years old
#' USA1990 = USA[USA$Year == 1990,]
#' x = 0:90
#' Ex = USA1990$Ex.Male[x+1]
#' Dx = USA1990$Dx.Male[x+1]
#' y = log(Dx/Ex)
#'
#'
#' ## Fit poisson and lognormal HP model and DLM
#' fit = hp(x = x, Ex = Ex, Dx = Dx, model = "poisson",
#'          M = 2000, bn = 1000, thin = 1)
#' fit2 = dlm(y, M = 100)
#'
#' ## Plot multiples life tables
#' plot(list(fit, fit2),
#'      age = 0:110, Ex = USA1990$Ex.Male,
#'      plotIC = FALSE, colors = c("red", "blue"),
#'      labels = c("HP Poisson", "DLM"))
#'
#'
#' ## Plot ClosedHP and ClosedDLM objects
#' close1 = hp_close(fit, method = "hp", x0 = 90)
#' close2 = dlm_close(fit2, method = "plateau")
#' plot(list(fit, fit2, close1, close2),
#'      plotIC = FALSE, colors = c("red", "blue", "green", "purple"),
#'      labels = c("HP", "DLM", "ClosedHP", "ClosedDLM"))
#'
#' @include fitted_dlm.R
#' @include fitted_hp.R
#' @include fun_aux.R
#'
#' @import ggplot2
#' @import scales
#'
#' @seealso [plot.DLM()], [plot.HP()], [plot.BLC()] and [plot.PredBLC()] for single plots.
#'
#' @export
plot.list <- function(x, age = NULL, Ex = NULL, plotIC = TRUE,
                      plotData = TRUE, labels = NULL, colors = NULL,
                      linetype = NULL, prob = 0.95, ...){
  fits = x
  if(length(fits) == 1){ stop("The length of fits must be two or more.") }
  if(all(unlist(lapply(fits, class)) %in% c("DLM", "ClosedDLM", "HP", "ClosedHP"))){

    ## Number of fitted models
    n_models <- length(fits)
    classes <- unlist(lapply(fits, class))
    h = rep(0, n_models)

    ## Checking ages
    if(is.null(age)) {
      age <- list(); data_aux <- list()
      for(i in 1:n_models){
        if(classes[i] %in% c("HP", "ClosedHP")){
          age[[i]] = fits[[i]]$data$x; data_aux[[i]] = fits[[i]]$data
        }else{
          age[[i]] = fits[[i]]$info$ages; data_aux[[i]] = data.frame(x = fits[[i]]$info$ages, qx = 1 - exp(-exp(fits[[i]]$info$y)))
        }
      }
    }else{
      warn = F; aux2 = age; ages_new = aux = age = list(); data_aux <- list()
      for(i in 1:n_models){
        if(classes[i] %in% c("HP", "ClosedHP")){
          age[[i]] = aux2
          data_aux[[i]] = fits[[i]]$data[fits[[i]]$data$x %in% age[[i]], ]
        }else if(classes[i] == "ClosedDLM"){
          if(any(!(aux2 %in% fits[[i]]$info$ages))) { warn = T }
          ages_new[[i]] = aux2[which(aux2 %in% fits[[i]]$info$ages)]
          age[[i]] = ages_new[[i]]

          ## selecting just the columns of the ages specified by the user
          aux[[i]] = which(fits[[i]]$info$ages %in% ages_new[[i]])

          fits[[i]]$qx = fits[[i]]$qx[,aux[[i]]]

          fits[[i]]$info$y = c(fits[[i]]$info$y[aux[[i]]])
          fits[[i]]$info$ages = fits[[i]]$info$ages[aux[[i]]]
          data_aux[[i]] = data.frame(x = fits[[i]]$info$ages, qx = 1-exp(-exp(fits[[i]]$info$y)))
        }else{
          if(min(aux2) < min(fits[[i]]$info$ages)) { warn = T }
          ages_new[[i]] = aux2[aux2 >= min(fits[[i]]$info$ages)]
          ages_to_predict = aux2[aux2 > max(fits[[i]]$info$ages)]; h[i] = length(ages_to_predict)
          age[[i]] = ages_new[[i]]

          ## selecting just the columns of the ages specified by the user
          aux[[i]] = which(fits[[i]]$info$ages %in% ages_new[[i]])

          ## The update is different according to object class
          fits[[i]]$mu = fits[[i]]$mu[,aux[[i]]]
          fits[[i]]$beta = fits[[i]]$beta[,aux[[i]]]

          fits[[i]]$info$y = c(fits[[i]]$info$y[aux[[i]]])
          fits[[i]]$info$ages = fits[[i]]$info$ages[aux[[i]]]
          data_aux[[i]] = data.frame(x = fits[[i]]$info$ages, qx = 1-exp(-exp(fits[[i]]$info$y)))
        }
      }
      if(warn){ warning("There are ages especified smaller than the ones in DLM fitted. These ages will not be used.") }
      rm(ages_new)
    }

    ####################################################################################

    ## qx fit and ci
    qx_fit = qx_cin = n_aux = NULL
    for(i in 1:n_models){
      aux = fitted(fits[[i]], age = age[[i]], Ex = Ex, prob = prob)
      if(h[i] > 0) {
        qx_pred <- predict(fits[[i]], h = h[i])
        aux_last_age = max(aux$age)
        aux_qx_fit = c(aux$qx.fitted, qx_pred$qx.fitted)
        aux_qi_fit = c(aux$qx.lower, qx_pred$qx.lower)
        aux_qs_fit = c(aux$qx.upper, qx_pred$qx.upper)
        aux1 = data.frame(age = c(aux$age, (aux_last_age+1):(aux_last_age+h[i])),
                          qx.fitted = aux_qx_fit)
        aux2 = data.frame(age = c(aux$age, (aux_last_age+1):(aux_last_age+h[i])),
                          qx.lower = aux_qi_fit, qx.upper = aux_qs_fit)
      }else{
        aux1 = data.frame(age = aux$age, qx.fitted = aux$qx.fitted)
        aux2 = data.frame(age = aux$age, qx.lower = aux$qx.lower, qx.upper = aux$qx.upper)
      }
      qx_fit <- rbind(qx_fit, aux1)
      qx_cin <- rbind(qx_cin, aux2)
      n_aux[i] <- nrow(aux1)
    }
    qx_fit$Model = paste("Model", rep(1:n_models, n_aux))
    qx_cin$Model = paste("Model", rep(1:n_models, n_aux))
    qx_fit = na.omit(qx_fit)
    qx_cin = na.omit(qx_cin)

    oneModel = F

    ## Customizing the plot
    if(n_models == 1){
      oneModel = T
      if(is.null(colors)) { colors = "seagreen" }
      if(is.null(labels)) { labels = "fitted" }
      if(is.null(linetype)) { linetype = "solid" }
    }

    ## checking if there are color or label inputs
    if(is.null(labels)){ labels <- unique(qx_fit$Model) }

    if(length(labels) != n_models) {
      warning("The number of labels does not match the number of models to plot.")
      labels <- unique(qx_fit$Model)
    }

    if(is.null(colors)){ colors = scales::hue_pal()(n_models) }

    if(length(colors) != n_models) {
      warning("The number of selected colors does not match the number of models to plot.")
      colors = scales::hue_pal()(n_models)
    }

    if(is.null(linetype)){ linetype = "solid" }
    if(length(linetype) == 1){ linetype = rep(linetype, n_models) }
    if(length(linetype) != n_models) {
      warning("The number of selected linetype must be one type or match the number of models to plot.")
      linetype = rep("solid", n_models)
    }

    ## Organizing data
    if(!oneModel){
      data = NULL
      for(i in 1:n_models){
        data_aux[[i]]$Model = paste("Model", i)
        data <- rbind(data, data_aux[[i]][,c("x", "qx", "Model")])
      }
      new_labels = labels
      new_colors = colors
      data <- na.omit(data)
      data$Model <- factor(data$Model, labels = new_labels, levels = paste("Model", 1:n_models))
      qx_cin$Model <- factor(qx_cin$Model, labels = new_labels, levels = paste("Model", 1:n_models))
      qx_fit$Model <- factor(qx_fit$Model, labels = new_labels, levels = paste("Model", 1:n_models))
    }else{
      data = data_aux[[1]][,c("x", "qx", "Model")]
      data$Model = "data"
      data <- na.omit(data)
      if(plotData){
        new_labels <- append(labels, "data", 0)
        new_colors <- append(colors, "gray10", 0)
      }else{
        new_labels = labels
        new_colors = colors
      }
    }
    ## lower limit:
    limits_y <- decimal(min(c(qx_cin$qx.lower[qx_cin$qx.lower > 0], data$qx[data$qx > 0], qx_fit$qx.fitted[qx_fit$qx.fitted > 0]), na.rm = T))

    ## Plot base:
    g <- ggplot2::ggplot() +
      ggplot2::scale_y_continuous(trans = 'log10', breaks = 10^-seq(limits_y,0),
                                  limits = 10^-c(limits_y,0), labels = scales::comma) +
      ggplot2::scale_x_continuous(breaks = seq(0, 200, by = 10), limits = c(NA, NA)) +
      ggplot2::xlab("Age") + ggplot2::ylab("qx") + ggplot2::theme_bw() +
      ggplot2::theme(plot.title = ggplot2::element_text(lineheight = 1.2),
                     axis.title.x = ggplot2::element_text(color = 'black', size = 12),
                     axis.title.y = ggplot2::element_text(color = 'black', size = 12),
                     axis.text = ggplot2::element_text(color = 'black', size = 12),
                     legend.text = ggplot2::element_text(size = 12),
                     legend.position = "bottom")

    if(plotIC){
      if(plotData){
        g + ggplot2::geom_point(data = data, ggplot2::aes(x = x, y = qx, col = Model), alpha = 0.8, size = 0.8) +
          ggplot2::geom_ribbon(data = qx_cin, ggplot2::aes(x = age, ymin = qx.lower, ymax = qx.upper, fill = Model), alpha = 0.3) +
          ggplot2::geom_line(data = qx_fit, ggplot2::aes(x = age, y = qx.fitted, col = Model, lty = Model), linewidth = 0.8, alpha = 0.8) +
          ggplot2::scale_colour_manual(name = NULL, values = new_colors, labels = new_labels) +
          ggplot2::scale_fill_manual(name = NULL, values = colors) +
          ggplot2::scale_linetype_manual(name = NULL, values = linetype, labels = new_labels) +
          ggplot2::guides(fill = "none")
      }else{
        g +
          ggplot2::geom_ribbon(data = qx_cin, ggplot2::aes(x = age, ymin = qx.lower, ymax = qx.upper, fill = Model), alpha = 0.3) +
          ggplot2::geom_line(data = qx_fit, ggplot2::aes(x = age, y = qx.fitted, col = Model, lty = Model), linewidth = 0.8, alpha = 0.8) +
          ggplot2::scale_colour_manual(name = NULL, values = new_colors, labels = new_labels) +
          ggplot2::scale_fill_manual(name = NULL, values = colors) +
          ggplot2::scale_linetype_manual(name = NULL, values = linetype, labels = new_labels) +
          ggplot2::guides(fill = "none")
      }
    }else{
      if(plotData){
        g + ggplot2::geom_point(data = data, ggplot2::aes(x = x, y = qx, col = Model), alpha = 0.8, size = 0.8) +
          ggplot2::geom_line(data = qx_fit, ggplot2::aes(x = age, y = qx.fitted, col = Model, lty = Model), linewidth = 0.8, alpha = 0.8) +
          ggplot2::scale_colour_manual(name = NULL, values = new_colors, labels = new_labels) +
          ggplot2::scale_linetype_manual(name = NULL, values = linetype, labels = new_labels) +
          ggplot2::guides(fill = "none")
      }else{
        g +
          ggplot2::geom_line(data = qx_fit, ggplot2::aes(x = age, y = qx.fitted, col = Model, lty = Model), linewidth = 0.8, alpha = 0.8) +
          ggplot2::scale_colour_manual(name = NULL, values = new_colors, labels = new_labels) +
          ggplot2::scale_linetype_manual(name = NULL, values = linetype, labels = new_labels) +
          ggplot2::guides(fill = "none")
      }
    }
  }else{
    stop("fits argument must be an object or a list of DLM and/or HP objects.")
  }
}
