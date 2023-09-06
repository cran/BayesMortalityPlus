#' @name improvement
#' @rdname improvement
#'
#' @title BLC: Improvement
#'
#' @description Calculates the improvement of each age, based on the resulting chains of the beta parameter from a fitted blc model.
#'
#' @usage
#' improvement(obj, prob = 0.95)
#'
#' @param obj A `BLC` object, result of a call to blc() function.
#' @param prob A real number that represents the credibility level of the intervals.
#'
#' @return A data.frame with the improvement values of each age, as well as their credible intervals.
#'
#' @examples
#' ## Importing log-mortality data from Portugal:
#' data(PT)
#' Y <- PT
#'
#' ## Fitting the model
#' fit = blc(Y = Y, M = 100, bn = 20)
#'
#' ## Improvement:
#' improvement(fit)
#' improvement(fit, prob = 0.9) #90% credible intervals
#'
#'
#' @export
improvement <- function(obj, prob = 0.95) {
	objClass <- class(obj)
	supportedClasses <- c("BLC", "ARBLC")

	if (!any(objClass %in% supportedClasses)) {
		stop("Invalid object type")
	}

	lower <- (1-prob)/2
	upper <- (1+prob)/2

	improvement <- apply(obj$beta[ ,(obj$bn+1):obj$M], 1, function (x) mean(1 - exp(-x)))
	improvement.int <- apply(obj$beta[ ,(obj$bn+1):obj$M], 1, function (x) quantile(1 - exp(-x),
	                                                                                        probs = c(lower, upper)))

	data.frame(improvement = improvement, lower.lim = improvement.int[1,], upper.lim = improvement.int[2,])
}
