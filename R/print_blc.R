#' @name print.BLC
#' @rdname print.BLC
#'
#' @title Print Values for BLC fitted models
#'
#' @description Print details from a fitted BLC model and returns it invisibly.
#'
#'
#' @param x A `BLC` object, result of a call to blc() function.
#' @param ... Further arguments passed to or from other methods.
#'
#' @return A character vector with the details of a fitted `BLC` model.
#'
#'
#' @seealso [print.DLM()], [print.HP()] and [print.PredBLC()] for `DLM`, `HP` or `PredBLC` methods.
#'
#' @export
print.BLC <- function(x, ...) {
  obj = x
	catf <- function(x, ...) cat(sprintf(x, ...), end = '\n')

	catf("Bayesian Lee-Carter Estimation\n")
	catf("Age groups: %d", nrow(obj$Y))
	catf("Time length: %d", ncol(obj$Y))
	catf("Iterations: %d (%d for warm-up)", obj$numit, obj$warmup)
	catf("Prior: N(%.2f, %.2e)", obj$m0, obj$C0)

	invisible(obj)
}
