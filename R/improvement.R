#' @name improvement
#' @rdname improvement
#'
#' @title Improvement for BLC fitted models
#'
#' @description Calculates the improvement of each age, based on the resulting chains of the beta parameter from a fitted blc model.
#'
#' @usage
#' improvement(obj, cred = 0.95)
#'
#' @param obj A `BLC` object, result of a call to blc() function.
#' @param cred A real number that represents the credibility level of the intervals.
#'
#' @return A data.frame with the improvement values of each age, as well as their credible intervals.
#'
#' @examples
#' ## Importing log-mortality data from Portugal:
#' data(PT)
#' Y <- PT
#'
#' ## Fitting the model
#' fit = blc(Y = Y, numit = 100, warmup = 20)
#'
#' ## Improvement:
#' improvement(fit)
#' improvement(fit, cred = 0.9) #90% credible intervals
#'
#'
#' @export
improvement <- function(obj, cred = 0.95) {
	objClass <- class(obj)
	supportedClasses <- c("BLC", "ARBLC")

	if (!any(objClass %in% supportedClasses)) {
		stop("Invalid object type")
	}

	lower <- (1-cred)/2
	upper <- (1+cred)/2

	improvement <- apply(obj$beta[ ,(obj$warmup+1):obj$numit], 1, function (x) mean(1 - exp(-x)))
	improvement.int <- apply(obj$beta[ ,(obj$warmup+1):obj$numit], 1, function (x) quantile(1 - exp(-x),
	                                                                                        probs = c(lower, upper)))

	data.frame(improvement = improvement, lower.lim = improvement.int[1,], upper.lim = improvement.int[2,])
}
