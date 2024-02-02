#' @name print.BLC
#' @rdname print.BLC
#'
#' @title BLC: Print
#'
#' @description Print details from a fitted BLC model and returns it invisibly.
#'
#'
#' @param x A `BLC` or `PredBLC`object, result of a call to blc() or predict() function.
#' @param ... Further arguments passed to or from other methods.
#'
#' @return A character vector with the details of a fitted `BLC` or `PredBLC` model.
#'
#'
#' @seealso [print.DLM()] and [print.HP()] for `DLM` or `HP` methods.
#'
#' @export
print.BLC <- function(x, ...) {
  obj = x
	catf <- function(x, ...) cat(sprintf(x, ...), end = '\n')

	catf("Bayesian Lee-Carter Estimation\n")
	catf("Age groups: %d", nrow(obj$Y))
	catf("Time length: %d", ncol(obj$Y))
	catf("Sample size: %d", obj$M)
	catf("Prior: N(%.2f, %.2e)", obj$m0, obj$C0)

	invisible(obj)
}

#' @export
print.PredBLC <- function(x, ...) {
  obj = x
  cat(sprintf('Forecast of a Bayesian Lee-Carter model (h = %d)\n',
              obj$h))
  cat("\n")

  invisible(obj)
}
