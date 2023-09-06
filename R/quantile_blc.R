#' @name quantile.BLC
#' @rdname quantile.BLC
#'
#' @title BLC: Sample quantiles
#'
#' @description Compute the quantiles based on the resulting chains from a fitted BLC model.
#'
#'
#' @param x A `BLC` object, result of a call to blc() function.
#' @param q A real number that represents the probability of the quantiles.
#' @param name A character with a parameter name of the blc model that should be returned. It can be one of these: "alpha", "beta", "kappa", "phiv", "theta", "phiw".
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
#' ## Parameters' median and quantiles 0.05, 0.95
#' quantile(fit, c(0.05, 0.5, 0.95), "alpha")
#' quantile(fit, c(0.05, 0.5, 0.95), "beta")
#' quantile(fit, c(0.05, 0.5, 0.95), "kappa")
#' quantile(fit, c(0.05, 0.5, 0.95), "phiv") ## random error precision
#' quantile(fit, c(0.05, 0.5, 0.95), "theta") ## drift parameter
#' quantile(fit, c(0.05, 0.5, 0.95), "phiw")
#'
#' @seealso [quantile.PredBLC()] for `PredBLC` method.
#'
#' @export
quantile.BLC <- function(x, q, name, ...) {
  obj = x
	matrixNames <- c("alpha", "beta", "kappa", "phiv")
	vectorNames <- c("theta", "phiw")

	if (name %in% matrixNames) {
		apply(obj[[name]][ ,(obj$bn+1):obj$M], 1, quantile, q)
	} else if (name %in% vectorNames) {
		quantile(obj[[name]][(obj$bn+1):obj$M], q)
	} else {
		stop("Invalid `name` argument")
	}
}
