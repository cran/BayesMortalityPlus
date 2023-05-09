#' @name fitted.DLM
#' @rdname fitted.DLM
#'
#' @title DLM: Fitted death probabilities (qx)
#'
#' @description This function computes the point estimations of the death probabilities (qx) of a
#' mortality graduation returned by dlm() or dlm_close() functions.
#'
#'
#' @param object Object of the following classes: `DLM` or `ClosedDLM`.
#' @param age Vector with the ages to calculate the death probabilities (Optional). By default, all ages are considered.
#' @param ... Other arguments.
#'
#' @return A data.frame object with the selected ages and the corresponding death probabilities.
#'
#' @examples
#' ## Importing mortality data from the USA available on the Human Mortality Database (HMD):
#' data(USA)
#'
#' ## Selecting the log mortality rate of the year 2000, ranging from 0 to 100 years old:
#' USA2000 = USA[USA$Year == 2000,]
#' x = 0:100
#' Ex = USA2000$Ex.Total[x+1]
#' Dx = USA2000$Dx.Total[x+1]
#'
#' qx_t = Dx/Ex
#' qx_t = 1 - exp(-qx_t)
#' y = log(qx_t)
#'
#' ## Fitting dlm
#' fit = dlm(y, M = 100, bn = 20, thin = 1)
#'
#' ## Estimating the death probabilities (qx)
#' fitted(fit)
#'
#' @seealso [fitted.HP()] and [fitted.BLC()] for `HP` or `BLC` methods.
#'
#' @export
fitted.DLM <- function(object, age = NULL, ...){
  fit = object
  qx_fitted = exp(fit$mu)
  qx_fitted = apply(qx_fitted, 2, median, na.rm = T)
  qx_fitted[(qx_fitted < 0 | qx_fitted > 1)] = NA
  df = data.frame(age = fit$info$ages, qx_fitted = qx_fitted)
  if(!is.null(age)) df = df[(df$age %in% age), ]
  return(df[order(df$age), ])
}

#' @export
fitted.ClosedDLM <- function(object, age = NULL, ...){

  fit = object
  qx_fitted = fit$qx
  qx_fitted[(qx_fitted < 0 | qx_fitted > 1)] = NA
  qx_fitted = apply(qx_fitted, 2, median, na.rm = T)

  close_age = fit$info$ages

  df = data.frame(age = close_age, qx_fitted = qx_fitted)
  if(!is.null(age)) df = df[(df$age %in% age), ]

  return(df[order(df$age), ])
}
