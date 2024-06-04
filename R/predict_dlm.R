#' @name predict.DLM
#' @rdname predict.DLM
#'
#' @title DLM: Prediction of death probability
#'
#' @description Extrapolates the mortality curve fitted by DLM by calculating the median
#' of death probability and the respective prediction interval.
#'
#'
#' @param object A `DLM` object that is result of a call to dlm() function.
#' @param h The ages prediction horizon.
#' @param prob Coverage probability of the predictive intervals.
#' @param ... Other arguments.
#'
#' @return A data.frame with the death probability prediction and prediction interval for the ages in the prediction horizon.
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
#' y = log(Dx/Ex)
#'
#' ## Fitting dlm
#' fit = dlm(y, M = 100)
#'
#' ## Extrapolating the death probabilities (qx)
#' predict(fit, h = 3, prob = 0.95)
#'
#'
#' @importFrom MASS mvrnorm
#'
#' @seealso [fitted.DLM()].
#'
#' @include ffbs.R
#'
#' @export
predict.DLM <- function(object, h, prob = 0.95, ...){

  fit = object
  N = length(fit$info$y)
  p = length(fit$info$prior$m0)
  y = fit$info$y
  Gt = fit$info$Gt
  Ft = fit$info$Ft
  delta = fit$info$delta[length(fit$info$delta)]
  # V = fit$sig2
  V = 0.01 ## same value used in Filtering
  sig2 = fit$sig2
  n = length(fit$sig2)

  aux = fit$param

  sim <- matrix(NA_real_, nrow = n, ncol = h)

  Wt = aux$Ct[N,,] * (1 - delta) / delta
  at = Gt %*% aux$mt[N,]
  Rt = Gt %*% aux$Ct[N,,] %*% t(Gt) + Wt
  ft = Ft %*% at
  Qt = (Ft %*% Rt %*% t(Ft) + V)[1,1]
  At = Rt %*% t(Ft) / Qt
  Ct = Rt - At %*% Ft %*% Rt   ## second moment

  # sim[, 1] <- rnorm(n, mean = ft, sd = sd(Qt * sig2))
  sim[, 1] <- rt(n, df = 2*aux$alpha)*sqrt(Qt*(aux$beta/aux$alpha)) + c(ft)

  if(h > 1) for(k in 2:h){
    #Wt = Ct
    at = Gt %*% at
    Rt = Gt %*% Rt %*% t(Gt) + Wt
    ft = Ft %*% at
    Qt = (Ft %*% Rt %*% t(Ft) + V)[1,1]
    At = Rt %*% t(Ft) / Qt
    Ct = Rt - At %*% Ft %*% Rt   ## second moment

    # sim[, k] <- rnorm(n, mean = ft, sd = sd(Qt * sig2))
    sim[, k] <- rt(n, df = 2*aux$alpha)*sqrt(Qt*(aux$beta/aux$alpha)) + c(ft)

  }

  qx_sim = exp(sim)
  qx_fitted = apply(qx_sim, 2, median, na.rm = T)
  qx_lim = apply(qx_sim, 2, quantile, probs = c((1-prob)/2, (1+prob)/2), na.rm = T)

  qx_fitted = data.frame(Ages = (fit$info$ages[N]+1):(fit$info$ages[N]+h), qx_fitted = qx_fitted)
  ret = data.frame(age = qx_fitted$Ages, qx.fitted = 1 - exp(-qx_fitted$qx_fitted),
                   qx.lower = 1 - exp(-qx_lim[1,]), qx.upper = 1 - exp(-qx_lim[2,]))
  ret[ret[,2:4] < 0,2:4] = 0
  ret[ret[,2:4] > 1,2:4] = 1

  return(ret)
}

