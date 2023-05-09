#' @name fitted.HP
#' @rdname fitted.HP
#'
#' @title HP Model: Fitted death probabilities (qx)
#'
#' @description This function computes the point estimations of the death probabilities (qx) of the `HP` or the `ClosedHP` class object fitted by the hp() or hp_close() functions.
#'
#'
#' @param object Object of the class `HP` or `ClosedHP` adjusted by the hp() or hp_close() functions.
#' @param age Vector with the ages to calculate the death probabilities (Optional). By default, all ages are considered.
#' @param ... Other arguments.
#'
#' @return A data.frame object with the selected ages and the corresponding death probabilities.
#'
#' @examples
#' ## Importing mortality data from the USA available on the Human Mortality Database (HMD):
#' data(USA)
#'
#' ## Selecting the exposure and the death count of the year 2000, ranging from 0 to 90 years old:
#' USA2000 = USA[USA$Year == 2000,]
#' x = 0:90
#' Ex = USA2000$Ex.Total[x+1]
#' Dx = USA2000$Dx.Total[x+1]
#'
#' ## Fitting a simple model:
#' fit = hp(x = x, Ex = Ex, Dx = Dx, M = 5000, bn = 0, thin = 10)
#'
#' ## Estimating the death probabilities (qx)
#' fitted(fit)
#' fitted(fit, age = 0:110)
#'
#' @include fun_aux.R
#'
#' @seealso [fitted.BLC()] and [fitted.DLM()] for `BLC` or `DLM` methods.
#'
#' @export
fitted.HP <- function(object, age = NULL, ...){
  fit = object
  if(is.null(age)){
    ## If the age argument is null, the ages in the fit model are used.
    age = fit$data$x
  }
  ## Aux object to save the qx markov chains
  qx_fitted = matrix(NA, nrow = nrow(fit$post.samples$mcmc_theta), ncol = length(age))

  if(fit$info$model %in% c("binomial","poisson")){
    ## point estimated qxs in each MCMC iteration
    for(i in 1:nrow(fit$post.samples$mcmc_theta)){
      qx_fitted[i,] = 1 - exp(-hp_curve_9(age, fit$post.samples$mcmc_theta[i,]))
    }
  }else{
    for(i in 1:nrow(fit$post.samples$mcmc_theta)){
      hp <- hp_curve(age, fit$post.samples$mcmc_theta[i,])
      qx_fitted[i,] = hp/(1 + hp)
    }
  }

  ## filtering invalid qx values (may happen in binomial models)
  qx_fitted[(qx_fitted < 0 | qx_fitted > 1)] = NA

  qx_fitted = apply(qx_fitted, 2, median, na.rm = T)
  return(data.frame(age = age, qx_fitted = qx_fitted))
}

#' @export
fitted.ClosedHP <- function(object, age = NULL, ...){

  fit = object
  qx_fitted = fit$qx
  qx_fitted[(qx_fitted < 0 | qx_fitted > 1)] = NA
  qx_fitted = apply(qx_fitted, 2, median)

  close_age = fit$data$x

  df = data.frame(age = close_age, qx_fitted = qx_fitted)

  if(!is.null(age)){

    df = df[(close_age %in% age), ]

    if(any(!(age %in% close_age))){
      age_not_fitted = age[!(age %in% close_age)]
      aux = data.frame(age = age_not_fitted, qx_fitted = NA)
      df = rbind(df, aux); row.names(df) = NULL
    }

  }
  return(df[order(df$age), ])
}
