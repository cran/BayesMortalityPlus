#' @name plot.HP
#' @rdname plot.HP
#'
#' @title HP: Plot the life table
#'
#' @description Function that returns a log-scale ggplot of `HP` and `ClosedHP` objects returned by the hp() and hp_close() functions.
#'
#'
#' @param x Object of the class `HP` or `ClosedHP` returned by hp() or hp_close() functions.
#' @param age Vector with the ages to plot the life table.
#' @param Ex Vector with the exposures of the selected ages. Its length must be equal to the age vector. This argument is only necessary when using poisson and binomial HP models.
#' @param plotIC Logical. If 'TRUE' (default), shows the predictive intervals.
#' @param plotData Logical. If 'TRUE' (default), shows crude rate (black dots).
#' @param labels Vector with the name of the curve label. (Optional).
#' @param colors Vector with the color of the curve. (Optional).
#' @param linetype Vector with the line type of the curve. (Optional).
#' @param prob Coverage probability of the predictive intervals. Default is '0.95'.
#' @param ... Other arguments.
#'
#' @return A 'ggplot' object with fitted life table.
#'
#' @examples
#' ## Selecting the exposure and the death count of the year 1990, ranging from 0 to 90 years old:
#' USA1990 = USA[USA$Year == 1990,]
#' x = 0:90
#' Ex = USA1990$Ex.Male[x+1]
#' Dx = USA1990$Dx.Male[x+1]
#'
#' ## Fitting the poisson and the lognormal model:
#' fit = hp(x = x, Ex = Ex, Dx = Dx, model = "poisson",
#'          M = 2000, bn = 1000, thin = 1)
#' fit2 = hp(x = x, Ex = Ex, Dx = Dx, model = "lognormal",
#'           M = 2000, bn = 1000, thin = 1)
#'
#' ## Plot the life tables:
#' plot(fit)
#' plot(fit2, age = 0:110, plotIC = TRUE)
#'
#' ## To plot multiples life tables see ?plot.list
#' plot(list(fit, fit2),
#'      age = 0:110, Ex = USA1990$Ex.Male,
#'      plotIC = FALSE, colors = c("red", "blue"),
#'      labels = c("Poisson", "Lognormal"))
#'
#' @include fitted_hp.R
#' @include fun_aux.R
#'
#' @import ggplot2
#' @import scales
#'
#' @seealso [plot.DLM()], [plot.BLC()] and [plot.PredBLC()] for `DLM`, `BLC` or `PredBLC` methods.
#'
#' [plot.list()] to the `list` method, adding multiple objects in one single plot.
#'
#' [plot_chain()] to plot the chains generated by the MCMC algorithms for the `HP` and `DLM` objects.
#'
#'
#' @export
plot.HP <- function(x, age = NULL, Ex = NULL, plotIC = TRUE,
                    plotData = TRUE, labels = NULL, colors = NULL,
                    linetype = NULL, prob = 0.95, ...){
  fit = x
  #### ----
  qx_fit = qx_ci = fitted(fit, age = age, Ex = Ex, prob = prob)

  ## Customizing the plot
  if(is.null(colors)) { colors = "seagreen" }
  if(is.null(labels)) { labels = "HP fitted" }
  if(is.null(linetype)) { linetype = "solid" }

  if(length(colors) != 1) {
    warning("colors length is incorrect. It will be replaced by default color.")
    colors = "seagreen"
  }
  if(length(labels) != 1) {
    warning("labels length is incorrect. It will be replaced by default label.")
    labels = "DLM fitted"
  }
  if(length(linetype) != 1) {
    warning("linetype length is incorrect. It will be replaced by default type.")
    labels = "solid"
  }

  ## Organizing data
  data = fit$data
  data$Model = "data"
  data <- na.omit(data)
  if(plotData){
    new_labels <- append(labels, "data", 0)
    new_colors <- append(colors, "gray10", 0)
    if(linetype %in% c("blank", "solid", "dashed", "dotted", "dotdash", "longdash", "twodash", "1F", "F1", "4C88C488", "12345678")){
      new_linetype = c("blank", linetype)
    }else{
      linetype = as.numeric(linetype)
      new_linetype = c(0, linetype)
    }
  }else{
    new_labels = labels
    new_colors = colors
  }

  ## lower limit:
  li_y <- decimal(min(c(qx_ci$qx.lower[qx_ci$qx.lower > 0], data$qx[data$qx > 0], qx_fit$qx.fitted[qx_fit$qx.fitted > 0]), na.rm = T))

  if(!is.null(age)) { data = data[data$x %in% age, ] }

  ## Plot base:
  g <- ggplot2::ggplot() +
    ggplot2::scale_y_continuous(trans = 'log10', breaks = 10^-seq(li_y,0), limits = 10^-c(li_y,0), labels = scales::comma) +
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
      g + ggplot2::geom_point(data = data, ggplot2::aes(x = x, y = qx, col = "data"), alpha = 1, size = 0.8) +
        ggplot2::geom_ribbon(data = qx_ci, ggplot2::aes(x = age, ymin = qx.lower, ymax = qx.upper, fill = "fitted"), alpha = 0.3) +
        ggplot2::geom_line(data = qx_fit, ggplot2::aes(x = age, y = qx.fitted, col = "fitted", lty = "fitted"), linewidth = 0.8, alpha = 0.8) +
        ggplot2::scale_colour_manual(name = NULL, values = new_colors, labels = new_labels) +
        ggplot2::scale_fill_manual(name = NULL, values = colors) +
        ggplot2::scale_linetype_manual(name = NULL, values = linetype) +
        ggplot2::guides(fill = "none", lty = "none",
                        color = ggplot2::guide_legend(override.aes = list(linetype = c(new_linetype),
                                                                          shape = c(19, NA))))
    }else{
      g +
        ggplot2::geom_ribbon(data = qx_ci, ggplot2::aes(x = age, ymin = qx.lower, ymax = qx.upper, fill = "fitted"), alpha = 0.3) +
        ggplot2::geom_line(data = qx_fit, ggplot2::aes(x = age, y = qx.fitted, col = "fitted", lty = "fitted"), linewidth = 0.8, alpha = 0.8) +
        ggplot2::scale_colour_manual(name = NULL, values = new_colors, labels = new_labels) +
        ggplot2::scale_fill_manual(name = NULL, values = colors) +
        ggplot2::scale_linetype_manual(name = NULL, values = linetype, labels = new_labels) +
        ggplot2::guides(fill = "none")
    }
  }else{
    if(plotData){
      g + ggplot2::geom_point(data = data, ggplot2::aes(x = x, y = qx, col = "data"), alpha = 1, size = 0.8) +
        ggplot2::geom_line(data = qx_fit, ggplot2::aes(x = age, y = qx.fitted, col = "fitted", lty = "fitted"), linewidth = 0.8, alpha = 0.8) +
        ggplot2::scale_colour_manual(name = NULL, values = new_colors, labels = new_labels) +
        ggplot2::scale_linetype_manual(name = NULL, values = linetype, labels = new_labels) +
        ggplot2::guides(fill = "none", lty = "none",
                        color = ggplot2::guide_legend(override.aes = list(linetype = c(new_linetype),
                                                                          shape = c(19, NA))))
    }else{
      g +
        ggplot2::geom_line(data = qx_fit, ggplot2::aes(x = age, y = qx.fitted, col = "fitted", lty = "fitted"), linewidth = 0.8, alpha = 0.8) +
        ggplot2::scale_colour_manual(name = NULL, values = new_colors, labels = new_labels) +
        ggplot2::scale_linetype_manual(name = NULL, values = linetype, labels = new_labels) +
        ggplot2::guides(fill = "none")
    }
  }
}

#' @export
plot.ClosedHP <- function(x, age = NULL, plotIC = TRUE, plotData = TRUE,
                          labels = NULL, colors = NULL, linetype = NULL,
                          prob = 0.95, ...){
  fit = x
  qx_fit = qx_ci = na.omit(fitted(fit, age = age, prob = prob))


  ## Customizing the plot
  ## Customizing the plot
  if(is.null(colors)) { colors = "seagreen" }
  if(is.null(labels)) { labels = "HP fitted" }
  if(is.null(linetype)) { linetype = "solid" }

  if(length(colors) != 1) {
    warning("colors length is incorrect. It will be replaced by default color.")
    colors = "seagreen"
  }
  if(length(labels) != 1) {
    warning("labels length is incorrect. It will be replaced by default label.")
    labels = "DLM fitted"
  }
  if(length(linetype) != 1) {
    warning("linetype length is incorrect. It will be replaced by default type.")
    labels = "solid"
  }

  ## Organizing data
  data = fit$data
  data$Model = "data"
  data <- na.omit(data)
  if(plotData){
    new_labels <- append(labels, "data", 0)
    new_colors <- append(colors, "gray10", 0)
    if(linetype %in% c("blank", "solid", "dashed", "dotted", "dotdash", "longdash", "twodash", "1F", "F1", "4C88C488", "12345678")){
      new_linetype = c("blank", linetype)
    }else{
      linetype = as.numeric(linetype)
      new_linetype = c(0, linetype)
    }
  }else{
    new_labels = labels
    new_colors = colors
  }

  ## lower limit:
  li_y <- decimal(min(c(qx_ci$qx.lower[qx_ci$qx.lower > 0], data$qx[data$qx > 0], qx_fit$qx.fitted[qx_fit$qx.fitted > 0]), na.rm = T))

  if(!is.null(age)) { data = data[data$x %in% age, ] }

  ## Plot base:
  g <- ggplot2::ggplot() +
    ggplot2::scale_y_continuous(trans = 'log10', breaks = 10^-seq(li_y,0), limits = 10^-c(li_y,0), labels = scales::comma) +
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
      g + ggplot2::geom_point(data = data, ggplot2::aes(x = x, y = qx, col = "data"), alpha = 1, size = 0.8) +
        ggplot2::geom_ribbon(data = qx_ci, ggplot2::aes(x = age, ymin = qx.lower, ymax = qx.upper, fill = "fitted"), alpha = 0.3) +
        ggplot2::geom_line(data = qx_fit, ggplot2::aes(x = age, y = qx.fitted, col = "fitted", lty = "fitted"), linewidth = 0.8, alpha = 0.8) +
        ggplot2::scale_colour_manual(name = NULL, values = new_colors, labels = new_labels) +
        ggplot2::scale_fill_manual(name = NULL, values = colors) +
        ggplot2::scale_linetype_manual(name = NULL, values = linetype) +
        ggplot2::guides(fill = "none", lty = "none",
                        color = ggplot2::guide_legend(override.aes = list(linetype = c(new_linetype),
                                                                          shape = c(19, NA))))
    }else{
      g +
        ggplot2::geom_ribbon(data = qx_ci, ggplot2::aes(x = age, ymin = qx.lower, ymax = qx.upper, fill = "fitted"), alpha = 0.3) +
        ggplot2::geom_line(data = qx_fit, ggplot2::aes(x = age, y = qx.fitted, col = "fitted", lty = "fitted"), linewidth = 0.8, alpha = 0.8) +
        ggplot2::scale_colour_manual(name = NULL, values = new_colors, labels = new_labels) +
        ggplot2::scale_fill_manual(name = NULL, values = colors) +
        ggplot2::scale_linetype_manual(name = NULL, values = linetype, labels = new_labels) +
        ggplot2::guides(fill = "none")
    }
  }else{
    if(plotData){
      g + ggplot2::geom_point(data = data, ggplot2::aes(x = x, y = qx, col = "data"), alpha = 1, size = 0.8) +
        ggplot2::geom_line(data = qx_fit, ggplot2::aes(x = age, y = qx.fitted, col = "fitted", lty = "fitted"), linewidth = 0.8, alpha = 0.8) +
        ggplot2::scale_colour_manual(name = NULL, values = new_colors, labels = new_labels) +
        ggplot2::scale_linetype_manual(name = NULL, values = linetype) +
        ggplot2::guides(fill = "none", lty = "none",
                        color = ggplot2::guide_legend(override.aes = list(linetype = c(new_linetype),
                                                                          shape = c(19, NA))))
    }else{
      g +
        ggplot2::geom_line(data = qx_fit, ggplot2::aes(x = age, y = qx.fitted, col = "fitted", lty = "fitted"), linewidth = 0.8, alpha = 0.8) +
        ggplot2::scale_colour_manual(name = NULL, values = new_colors, labels = new_labels) +
        ggplot2::scale_linetype_manual(name = NULL, values = linetype, labels = new_labels) +
        ggplot2::guides(fill = "none")
    }
  }
}
