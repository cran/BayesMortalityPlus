#' @name mean.BLC
#' @rdname mean.BLC
#'
#' @title BLC: Arithmetic mean
#'
#' @description Calculates the means based on the resulting chains from a fitted BLC model.
#'
#'
#' @param x A `BLC` object, result of a call to blc() function.
#' @param name A character with a parameter name of the BLC model that should be returned. It can be one of these: "alpha", "beta", "kappa", "phiv", "theta", "phiw".
#' @param ... Further arguments passed to or from other methods.
#'
#' @return A vector with the mean values of the selected parameter.
#'
#' @examples
#' data(PT)
#' Y <- PT
#'
#' ## Fitting the model
#' fit = blc(Y = Y, M = 100, bn = 20)
#'
#' mean(fit, "kappa")
#'
#' @seealso [mean.PredBLC()] for `PredBLC` object method.
#'
#'
#' @export
mean.BLC <- function(x, name, ...) {
  obj = x
	matrixNames <- c("alpha", "beta", "kappa", "phiv")
	vectorNames <- c("theta", "phiw")

	if (name %in% matrixNames) {
		apply(obj[[name]][ ,(obj$bn+1):obj$M], 1, mean)
	} else if (name %in% vectorNames) {
		mean(obj[[name]][(obj$bn+1):obj$M])
	} else {
		stop("Invalid `name` argument")
	}
}
