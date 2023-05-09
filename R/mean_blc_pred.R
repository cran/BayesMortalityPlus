#' @name mean.PredBLC
#' @rdname mean.PredBLC
#'
#' @title Arithmetic Mean for Predictions
#'
#' @description Calculates the means based on the resulting chains from a predicted year.
#'
#'
#' @param x A `PredBLC` object, result to the pred() function call on a `BLC` object.
#' @param h A positive integer specifying the year in the prediction horizon to be calculated.
#' @param ... Further arguments passed to or from other methods.
#'
#' @return A vector with the mean values of the log-mortality chains.
#'
#' @examples
#' data(PT)
#' Y <- PT
#'
#' ## Fitting the model
#' fit = blc(Y = Y, numit = 100, warmup = 20)
#'
#' ## Prediction for 2 years ahead
#' pred = predict(fit, h = 2)
#'
#' mean(pred, 1)
#' mean(pred, 2)
#'
#' @seealso [mean.BLC()] for `BLC` object method.
#'
#' @export
mean.PredBLC <- function(x, h, ...) {
  obj = x
	apply(obj$y[ ,h, ], 2, mean)
}
