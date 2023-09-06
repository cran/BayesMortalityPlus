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
#' mx_t = Dx/Ex
#' qx_t = 1 - exp(-mx_t)
#' y = log(qx_t)
#'
#' ## Fitting dlm
#' fit = dlm(y, M = 100, bn = 20, thin = 1)
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
  delta = fit$info$delta
  V = fit$sig2
  n = length(V)

  sim <- matrix(NA_real_, nrow = n, ncol = h)

  for(i in 1:n){
    aux = ff(m0 = fit$info$prior$m0, C0 = fit$info$prior$C0, y = as.matrix(y),
             V = V[i], Ft = Ft, Gt = Gt, delta = delta)

    Wt = aux$C[N,,] * (1 - delta) / delta
    at = Gt %*% aux$m[N,]
    Rt = Gt %*% aux$C[N,,] %*% t(Gt) + Wt
    ft = Ft %*% at
    Qt = Ft %*% Rt %*% t(Ft) + V[i]
    At = Rt %*% t(Ft) %*% chol2inv(chol(Qt))
    Ct = Rt - At %*% Ft %*% Rt   ## second moment

    sim[i, 1] <- MASS::mvrnorm(1, ft, Qt)

    if(h > 1) for(k in 2:h){
      Wt = Ct * (1 - delta) / delta
      at = Gt %*% at
      Rt = Gt %*% Rt %*% t(Gt) + Wt
      ft = Ft %*% at
      Qt = Ft %*% Rt %*% t(Ft) + V[i]
      At = Rt %*% t(Ft) %*% chol2inv(chol(Qt))
      Ct = Rt - At %*% Ft %*% Rt   ## second moment

      sim[i, k] <- MASS::mvrnorm(1, ft, Qt)
    }
  }

  qx_sim = exp(sim)
  qx_sim[qx_sim < 0] = 0
  qx_sim[qx_sim > 1] = 1
  qx_fitted = apply(qx_sim, 2, median, na.rm = T)
  qx_lim = apply(qx_sim, 2, quantile, probs = c((1-prob)/2, (1+prob)/2), na.rm = T)
  qx_fitted = data.frame(Ages = (fit$info$ages[N]+1):(fit$info$ages[N]+h), qx_fitted = qx_fitted)
  return(data.frame(Ages = qx_fitted$Ages, qx_fitted = qx_fitted$qx_fitted,
                    qx_inf = qx_lim[1,], qx_sup = qx_lim[2,]))
}

