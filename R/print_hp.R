#' @name print.HP
#' @rdname print.HP
#'
#' @title HP: Print
#'
#' @description Print details from a fitted `HP` or `ClosedHP` models and returns it invisibly.
#'
#'
#' @param x A `HP` or `ClosedHP` object, result of a call to hp() or hp_close() function.
#' @param ... Further arguments passed to or from other methods.
#'
#' @return A character vector with the details of a fitted `HP` or `ClosedHP` model.
#'
#' @examples
#'  ## Importing mortality data from the USA available on the Human Mortality Database (HMD):
#' data(USA)
#'
#' ## Selecting the exposure and death count of the 2010 male population ranging from 0 to 90 years old
#' USA2010 = USA[USA$Year == 2010,]
#' x = 0:90
#' Ex = USA2010$Ex.Male[x+1]
#' Dx = USA2010$Dx.Male[x+1]
#'
#' ## Fitting binomial model
#' fit = hp(x = x, Ex = Ex, Dx = Dx, M = 5000, bn = 0, thin = 10)
#' print(fit)
#'
#' @seealso [print.DLM()] and [print.BLC()] for `DLM` or `BLC` methods.
#'
#' @export
print.HP <- function(x, ...){
  fit = x
  catf <- function(x, ...) cat(sprintf(x, ...), end = '\n')

  catf("HP curve estimation")
  catf("Model: %s", fit$info$model)
}

#' @export
print.ClosedHP <- function(x, ...){
  fit = x
  catf <- function(x, ...) cat(sprintf(x, ...), end = '\n')

  catf("HP closure curve estimation")
  catf("Method: %s", fit$method)
  catf("Min. age: %s", min(fit$data$x, na.rm = T))
  catf("Max. age: %s", max(fit$data$x, na.rm = T))
}
