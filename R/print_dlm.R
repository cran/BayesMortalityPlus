#' @name print.DLM
#' @rdname print.DLM
#'
#' @title DLM: Print
#'
#' @description Print details from a fitted `DLM` or `ClosedDLM` models and returns it invisibly.
#'
#'
#' @param x A `DLM` or `ClosedDLM` object, result of a call to dlm() or dlm_close() function.
#' @param ... Further arguments passed to or from other methods.
#'
#' @return A character vector with the details of a fitted `DLM` or `ClosedDLM` model.
#'
#' @examples
#' ## Importing mortality data from the USA available on the Human Mortality Database (HMD):
#' data(USA)
#'
#' ## Selecting the log mortality rate of the 2010 male population ranging from 0 to 100 years old
#' USA2010 = USA[USA$Year == 2010,]
#' x = 0:100
#' Ex = USA2010$Ex.Male[x+1]
#' Dx = USA2010$Dx.Male[x+1]
#' y = log(Dx/Ex)
#'
#' ## Fitting DLM
#' fit = dlm(y, M = 100, bn = 20, thin = 1)
#' print(fit)
#'
#' @seealso [print.HP()] and [print.BLC()] for `HP` or `BLC` methods.
#'
#' @export
print.DLM <- function(x, ...){
  fit = x
  catf <- function(x, ...) cat(sprintf(x, ...), end = '\n')

  catf("DLM for life tables fitted")
  catf("")
  catf("Ft:")
  cat("[", fit$info$Ft,"]")
  catf("")
  catf("")
  catf("Gt:")
  for(i in 1:nrow(fit$info$Gt)) {cat(fit$info$Gt[i,]); catf("")}
  catf("")
  catf("Discount factor: %s", as.character(fit$info$delta))
  catf("")
  catf("Ages fitted:")
  cat(fit$info$ages)
  catf("")
}

#' @export
print.ClosedDLM <- function(x, ...){
  fit = x
  catf <- function(x, ...) cat(sprintf(x, ...), end = '\n')

  catf("DLM closure curve estimation")
  catf("Method: %s", fit$method)
  catf("Min. age: %s", min(fit$info$ages, na.rm = T))
  catf("Max. age: %s", max(fit$info$ages, na.rm = T))
}
