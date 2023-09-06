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
#' @param prob Coverage probability of the predictive intervals.
#' @param ... Other arguments.
#'
#' @return A data.frame object with the selected ages and the corresponding estimates and predictive intervals of the death probabilities.
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
#' y = log(Dx/Ex)
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
fitted.DLM <- function(object, age = NULL, prob = 0.95, ...){

  fit = object

  qx_fitted = 1 - exp(-exp(fit$mu))
  qx_fitted = apply(qx_fitted, 2, median, na.rm = T)
  qx_fitted[(qx_fitted < 0 | qx_fitted > 1)] = NA_real_

  t = ncol(fit$mu)
  n = nrow(fit$mu)
  fitted = matrix(NA_real_, nrow = n, ncol = t)

  # for(i in 1:n){
  #   sim = rnorm(t, fit$mu[i,], sqrt(fit$sig2[i]))
  #   fitted[i,] <- exp(sim)
  # }
  for(i in 1:t){
    sim = rnorm(n, fit$mu[,i], sqrt(fit$sig2))
    fitted[,i] <- 1 - exp(-exp(sim))
  }

  qi = apply(fitted, 2, quantile, (1-prob)/2, na.rm = T)
  qs = apply(fitted, 2, quantile, (1+prob)/2, na.rm = T)

  aux = data.frame(age = fit$info$ages, qx_fitted = qx_fitted, qi = qi, qs = qs)
  aux[!(aux$qi > 0), 3] = 0
  aux[!(aux$qs < 1), 4] = 1

  if(!is.null(age)) aux = aux[(aux$age %in% age), ]

  return(aux)
}

#' @export
fitted.ClosedDLM <- function(object, age = NULL, prob = 0.95, ...){

  fit = object

  fitted = fit$qx
  close_age = fit$info$ages

  qx_fitted = fitted
  qx_fitted[(qx_fitted < 0 | qx_fitted > 1)] = NA_real_
  qx_fitted = apply(qx_fitted, 2, median, na.rm = T)

  qi = apply(fitted, 2, quantile, (1-prob)/2, na.rm = T)
  qs = apply(fitted, 2, quantile, (1+prob)/2, na.rm = T)

  df = data.frame(age = close_age, qx_fitted = qx_fitted, qi = qi, qs = qs)
  df[!(df$qi > 0), 3] = 0
  df[!(df$qs < 1), 4] = 1

  if(!is.null(age)) df = df[(df$age %in% age), ]
  return(df[order(df$age), ])
}

