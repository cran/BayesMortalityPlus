#' @name expectancy
#' @rdname expectancy
#' @title Generic expectancy function
#'
#' @description Generic function to \code{expectancy} method.
#'
#'
#' @param x Object of one of these class: `HP`, `DLM`, `BLC`, `ClosedHP`, `ClosedDLM`, `BLC`, or `PredBLC`.
#' @param ... Further arguments passed to or from other methods.
#'
#'
#' @return
#' A data.frame and (if graph = TRUE) a plot for `HP`, `DLM`, `ClosedHP` and `ClosedDLM` methods.
#' A list that contains three vectors with the fitted values of life expectancy and the lower and upper limits of the credible intervals for each year used in fitted model or for the prediction, for `BLC` and `PredBLC` methods.
#'
#' @details  This function computes the life expectancy given by:
#'
#'  \eqn{e_x = \sum tp_x}
#'
#' where:
#'
#'  \eqn{tp_x = p_0 x p_1 x ... x p_x}
#'
#' @seealso  [expectancy.HP()], [expectancy.DLM()] and [expectancy.BLC()].
#' @export
expectancy = function(x, ...) UseMethod("expectancy")

#'
#' @export
expectancy.default = function(x, ...) print(x, ...)
