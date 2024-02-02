#' @name plot.PredBLC
#' @rdname plot.PredBLC
#'
#' @title BLC: Plot the log-mortality of a prediction
#'
#' @description This functions plot the predicted log-mortality and the predict
#' intervals of the log-mortality for a specific year in the prediction horizon
#'
#'
#' @param x A `PredBLC` object, result to the pred() function call on a `BLC` object.
#' @param h A numeric vector specifying the year(s) in the prediction horizon to be calculated.
#' @param prob A real number that represents the probability of the predict interval.
#' @param plotIC Logical. If 'TRUE' (default), shows the predictive intervals.
#' @param age A numeric vector indicating the modelled ages. (Optional).
#' @param ... Other arguments.
#'
#' @return A 'ggplot' object with the predicted mortality rates and their predict intervals.
#'
#' @examples
#' ## Importing log-mortality data from Portugal:
#' data(PT)
#' Y <- PT
#'
#' ## Fitting the model
#' fit = blc(Y = Y, M = 100, bn = 20)
#'
#' #' ## Prediction for 10 years ahead
#' pred = predict(fit, h = 3)
#'
#' ## Plotting
#' plot(pred, h = 1)
#' \donttest{plot(pred, h = 3, prob = 0.9)}
#'
#' @seealso [plot.HP()], [plot.DLM()] and [plot.BLC] for `HP`, `DLM` or `BLC` methods.
#'
#' @import ggplot2
#' @import scales
#'
#' @export
plot.PredBLC <- function(x, h = NULL, prob = 0.95, plotIC = TRUE, age = NULL,
                         ...) {

  obj = x
  alpha <- 1 - prob

  if(is.null(h))	{ h <- dim(obj$y)[2] }
  if(any(h > dim(obj$y)[2]))	{ stop("h has invalid values.") }

  h.size <- length(h)
  res <- array(dim = c(h.size, 3, dim(obj$y)[3]))

  for(ind in 1:h.size){
    for (i in 1:dim(obj$y)[3]) {
      tmp <- obj$y[ ,h[ind], i]
      res[ind, 1, i] <- mean(tmp)
      res[ind, 2, i] <- quantile(tmp, 1 - alpha/2)
      res[ind, 3, i] <- quantile(tmp, alpha/2)
    }
  }
  res <- exp(res)
  n <- length(res[1,1,])
  if(is.null(age)) {age = 1:n}

  if(h.size == 1){
    if(plotIC){
      ggplot2::ggplot(data = NULL) +
        ggplot2::scale_y_continuous(trans = 'log10', breaks = 10^-seq(20,0),
                                    limits = 10^-c(NA_real_,0), labels = scales::comma) +
        ggplot2::scale_x_continuous(breaks = seq(0, 200, by = 10), limits = c(NA_real_, NA_real_)) +
        ggplot2::xlab("x (index)") + ggplot2::ylab("qx") + ggplot2::theme_bw() +
        ggplot2::theme(plot.title = ggplot2::element_text(lineheight = 1.2),
                       axis.title.x = ggplot2::element_text(color = 'black', size = 12),
                       axis.title.y = ggplot2::element_text(color = 'black', size = 12),
                       axis.text = ggplot2::element_text(color = 'black', size = 12),
                       legend.text = ggplot2::element_text(size = 12),
                       legend.position = "bottom") +
        ggplot2::geom_ribbon(ggplot2::aes(x = age, ymin = res[1,3,], ymax = res[1,2,]), fill = "steelblue4", alpha = 0.3) +
        ggplot2::geom_line(ggplot2::aes(x = age, y = res[1,1,]), col = "steelblue", linewidth = 0.8, alpha = 0.8)
    }else{
      ggplot2::ggplot(data = NULL) +
        ggplot2::scale_y_continuous(trans = 'log10', breaks = 10^-seq(20,0),
                                    limits = 10^-c(NA_real_,0), labels = scales::comma) +
        ggplot2::scale_x_continuous(breaks = seq(0, 200, by = 10), limits = c(NA_real_, NA_real_)) +
        ggplot2::xlab("x (index)") + ggplot2::ylab("qx") + ggplot2::theme_bw() +
        ggplot2::theme(plot.title = ggplot2::element_text(lineheight = 1.2),
                       axis.title.x = ggplot2::element_text(color = 'black', size = 12),
                       axis.title.y = ggplot2::element_text(color = 'black', size = 12),
                       axis.text = ggplot2::element_text(color = 'black', size = 12),
                       legend.text = ggplot2::element_text(size = 12),
                       legend.position = "bottom") +
        ggplot2::geom_line(ggplot2::aes(x = age, y = res[1,1,]), col = "steelblue", linewidth = 0.8, alpha = 0.8)
    }

  }else{

    aux_res <- cbind(t(res[1,,]), h[1])
    for(ind in 2:h.size){ aux_res <- rbind(aux_res, cbind(t(res[ind,,]), h[ind])) }
    df_res <- data.frame(aux_res)
    index = rep(age, h.size)

    if(plotIC){
      ggplot2::ggplot(data = df_res) +
        ggplot2::scale_y_continuous(trans = 'log10', breaks = 10^-seq(20,0),
                                    limits = 10^-c(NA_real_,0), labels = scales::comma) +
        ggplot2::scale_x_continuous(breaks = seq(0, 200, by = 10), limits = c(NA_real_, NA_real_)) +
        ggplot2::xlab("x (index)") + ggplot2::ylab("qx") + ggplot2::theme_bw() +
        ggplot2::theme(plot.title = ggplot2::element_text(lineheight = 1.2),
                       axis.title.x = ggplot2::element_text(color = 'black', size = 12),
                       axis.title.y = ggplot2::element_text(color = 'black', size = 12),
                       axis.text = ggplot2::element_text(color = 'black', size = 12),
                       legend.text = ggplot2::element_text(size = 12),
                       legend.position = "bottom") + ggplot2::labs(color = NULL, fill = NULL) +
        ggplot2::geom_ribbon(ggplot2::aes(x = index, ymin = X3, ymax = X2, fill = paste0("h = ", as.factor(X4))), alpha = 0.3) +
        ggplot2::geom_line(ggplot2::aes(x = index, y = X1, col = paste0("h = ", as.factor(X4))), linewidth = 0.8, alpha = 0.8)
    }else{
      ggplot2::ggplot(data = df_res) +
        ggplot2::scale_y_continuous(trans = 'log10', breaks = 10^-seq(20,0),
                                    limits = 10^-c(NA_real_,0), labels = scales::comma) +
        ggplot2::scale_x_continuous(breaks = seq(0, 200, by = 10), limits = c(NA_real_, NA_real_)) +
        ggplot2::xlab("x (index)") + ggplot2::ylab("qx") + ggplot2::theme_bw() +
        ggplot2::theme(plot.title = ggplot2::element_text(lineheight = 1.2),
                       axis.title.x = ggplot2::element_text(color = 'black', size = 12),
                       axis.title.y = ggplot2::element_text(color = 'black', size = 12),
                       axis.text = ggplot2::element_text(color = 'black', size = 12),
                       legend.text = ggplot2::element_text(size = 12),
                       legend.position = "bottom") + ggplot2::labs(color = NULL, fill = NULL) +
        ggplot2::geom_line(ggplot2::aes(x = index, y = X1, col = paste0("h = ", as.factor(X4))), linewidth = 0.8, alpha = 0.8)
    }
  }
}

