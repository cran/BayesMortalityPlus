#' @name print.PredBLC
#' @rdname print.PredBLC
#'
#' @title Print Values for BLC prediction models
#'
#' @description Print details from a fitted BLC prediction model and returns it invisibly.
#'
#'
#' @param x A `PredBLC` object, result to the pred() function call on a `BLC` object.
#' @param ... Further arguments passed to or from other methods.
#'
#' @return A character vector with the details of a fitted `PredBLC` model.
#'
#'
#' @seealso [print.DLM()], [print.HP()] and [print.BLC()] for `DLM`, `HP` or `BLC` methods.
#'
#' @export
print.PredBLC <- function(x, ...) {
  obj = x
	cat(sprintf('Forecast of a Bayesian Lee-Carter model (h = %d)\n',
				obj$h))

	invisible(obj)
}
