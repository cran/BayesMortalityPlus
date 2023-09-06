#' @name quantile.PredBLC
#' @rdname quantile.PredBLC
#'
#' @title BLC: Sample quantiles for predictions
#'
#' @description Calculates the quantiles of log-mortality based on the resulting chains from a predicted year.
#'
#'
#' @param x A `PredBLC` object, result to the pred() function call on a `BLC` object.
#' @param q A real number that represents the probability of the quantiles.
#' @param h A positive integer especifying the year in the prediciton horizon to be calculated.
#' @param ... Further arguments passed to or from other methods.
#'
#' @return A data.frame with the quantiles of the selected parameter.
#'
#' @examples
#' ## Importing log-mortality data from Portugal:
#' data(PT)
#' Y <- PT
#'
#' ## Fitting the model
#' fit = blc(Y = Y, M = 100, bn = 20)
#'
#' ## Prediction for 2 years ahead
#' pred = predict(fit, h = 2)
#'
#' ## The log-mortality median for the first year of prediction
#' quantile(pred, q = 0.5, h = 1)
#'
#' ## The 0.1 and 0.9 quantiles for the first and second year of prediction
#' quantile(pred, q = c(0.1, 0.9), h = 1)
#' quantile(pred, q = c(0.1, 0.9), h = 2)
#'
#' @seealso [quantile.BLC()] for `BLC` method.
#'
#' @export
quantile.PredBLC <- function(x, q, h, ...) {
  obj = x
	if(length(h) > 1) h = h[1]
	apply(obj$y[ ,h, ], 2, quantile, q)
}
