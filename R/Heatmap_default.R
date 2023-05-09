#' @name Heatmap
#' @rdname Heatmap
#' @title Generic Heatmap function
#'
#' @description Generic function to \code{Heatmap} method.
#'
#' @param x Object or list of objects of class `HP`, `DLM`, `ClosedHP` or `ClosedDLM`. Object of class `BLC` or `PredBLC`.
#' @param ... Further arguments passed to or from other methods.
#'
#' @return A ggplot2 heatmap of the life expectancy.
#'
#' @seealso  [Heatmap.HP()], [Heatmap.DLM()], [Heatmap.BLC()] and [Heatmap.list()].
#' @export
Heatmap = function(x, ...) UseMethod("Heatmap")

#'
#' @export
Heatmap.default = function(x, ...) print(x, ...)
