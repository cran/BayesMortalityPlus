#' @name fitted.HP
#' @rdname fitted.HP
#'
#' @title HP: Fitted death probabilities (qx)
#'
#' @description This function computes the point estimations of the death probabilities (qx) of the `HP` or the `ClosedHP` class object fitted by the hp() or hp_close() functions.
#'
#'
#' @param object Object of the class `HP` or `ClosedHP` adjusted by the hp() or hp_close() functions.
#' @param age Vector with the ages to calculate the death probabilities (Optional). By default, all ages are considered.
#' @param Ex Vector with the exposures of the selected ages. Its length must be equal to the age vector. This argument is only necessary when using the Poisson and the Binomial distributions.
#' @param prob Coverage probability of the predictive intervals.
#' @param ... Other arguments.
#'
#' @return A data.frame object with the selected ages and the corresponding estimates and predictive intervals of the death probabilities.
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
#' fitted(fit, age = 0:110, Ex = USA2000$Ex.Total[0:110+1])
#'
#' @include fun_aux.R
#'
#' @seealso [fitted.BLC()] and [fitted.DLM()] for `BLC` or `DLM` methods.
#'
#' @export
fitted.HP <- function(object, age = NULL, Ex = NULL, prob = 0.95, ...){

  set.seed(123) ## Set seed to reproducibility
  fit = object

  ## checking if age and Ex were inputed by the user
  if(is.null(age) && is.null(Ex)){
    ## if age and Ex are null, fetch from the fit model
    age = fit$data$x
    Ex = fit$data$Ex

  }else if(is.null(age) && !(is.null(Ex))){

    if(fit$info$model %in% c("binomial", "poisson")) { stop("Missing age argument.") }
    age = fit$data$x

  }else if(!(is.null(age)) && is.null(Ex)){

    if(fit$info$model %in% c("binomial","poisson")){

      if(all(age %in% fit$data$x)){
        min_age = min(fit$data$x, na.rm = T)
        Ex = fit$data$Ex[age-min_age+1]
      }else{
        stop("Missing Ex argument.")
      }

    }
  }else if(length(age) != length(Ex)){
    ## length check for age and Ex
    stop("age and Ex arguments have different lengths.")
  }

  ## checking for invalid probabilities
  if(prob < 0 || prob > 1){ stop("Invalid death probability values.") }

  if(fit$info$model == "binomial"){

    age_out = age[is.na(Ex)]
    age = age[!is.na(Ex)] ## Removing ages with no exposures
    Ex = Ex[!is.na(Ex)] ## Removing exposures with NA values

    ## Aux object to save the qx markov chains
    qx_fitted = matrix(NA_real_, nrow = nrow(fit$post.samples$mcmc_theta), ncol = length(age))
    qx_ic = matrix(NA_real_, nrow = nrow(fit$post.samples$mcmc_theta), ncol = length(age))

    for (i in 1:nrow(qx_ic)){
      qx = 1 - exp(-hp_curve_9(age, fit$post.samples$mcmc_theta[i,]))
      qx = ifelse((qx < 0 | qx > 1), NA_real_, qx)
      qx_fitted[i,] = qx ## # Estimativa pontual
      sim = rbinom(length(age), trunc(Ex), qx)
      qx_ic[i,] = sim/Ex
    }

  }else if(fit$info$model == "poisson"){

    age_out = age[is.na(Ex)]
    age = age[!is.na(Ex)] ## Removing ages with no exposures
    Ex = Ex[!is.na(Ex)] ## Removing exposures with NA values

    ## Aux object to save the qx markov chains
    qx_fitted = matrix(NA_real_, nrow = nrow(fit$post.samples$mcmc_theta), ncol = length(age))
    qx_ic = matrix(NA_real_, nrow = nrow(fit$post.samples$mcmc_theta), ncol = length(age))

    for (i in 1:nrow(qx_ic)){
      qx = 1 - exp(-hp_curve_9(age, fit$post.samples$mcmc_theta[i,]))
      qx = ifelse((qx < 0 | qx > 1), NA_real_, qx)
      qx_fitted[i,] = qx ## # Estimativa pontual
      sim = rpois(length(age), lambda = Ex*qx)
      qx_ic[i,] = sim/Ex
    }

  }else{

    age_out = NULL
    ## Aux object to save the qx markov chains
    qx_fitted = matrix(NA_real_, nrow = nrow(fit$post.samples$mcmc_theta), ncol = length(age))
    qx_ic = matrix(NA_real_, nrow = nrow(fit$post.samples$mcmc_theta), ncol = length(age))

    for(i in 1:nrow(qx_ic)){
      hp <- hp_curve(age, fit$post.samples$mcmc_theta[i,])
      qx_fitted[i,] = hp/(1 + hp)
      sim = rnorm(length(age), log(hp), sqrt(fit$post.samples$sigma2[i]))
      qx_ic[i,] <- exp(sim)/(1+exp(sim))
    }

  }

  qi = apply(qx_ic, 2, quantile, (1-prob)/2, na.rm = T)
  qs = apply(qx_ic, 2, quantile, (1+prob)/2, na.rm = T)
  qx_fitted = apply(qx_fitted, 2, median, na.rm = T)

  aux = data.frame(age = age, qx.fitted = qx_fitted, qx.lower = qi, qx.upper = qs)
  aux[!(aux$qx.lower > 0), 3] = 0
  aux[!(aux$qx.upper < 1), 4] = 1

  if(length(age_out) > 0){
    aux2 <- data.frame(age = age_out, qx.fitted = NA_real_, qx.lower = NA_real_, qx.upper = NA_real_)
    aux <- rbind(aux, aux2)
    aux <- aux[order(aux$age),]
  }

  return(aux)
}

#' @export
fitted.ClosedHP <- function(object, age = NULL, prob = 0.95, ...){

  fit = object

  if(fit$method == "Mix"){
    qx_fitted = apply(fit$qx, 2, median)
    return(data.frame(age = fit$data$x, qx.fitted = qx_fitted, qx.lower = NA_real_, qx.upper = NA_real_))
  }

  qx_fitted = fit$qx
  qx_fitted[(qx_fitted < 0 | qx_fitted > 1)] = NA_real_
  close_age = fit$data$x

  qi = apply(qx_fitted, 2, quantile, (1-prob)/2, na.rm = T)
  qs = apply(qx_fitted, 2, quantile, (1+prob)/2, na.rm = T)
  qx_fitted = apply(qx_fitted, 2, median, na.rm = T)

  df = data.frame(age = close_age, qx.fitted = qx_fitted, qx.lower = qi, qx.upper = qs)
  df[!(df$qi > 0), 3] = 0
  df[!(df$qs < 1), 4] = 1

  if(!is.null(age)){

    df = df[(close_age %in% age), ]

    if(any(!(age %in% close_age))){
      age_not_fitted = age[!(age %in% close_age)]
      aux = data.frame(age = age_not_fitted, qx.fitted = NA_real_, qx.lower = NA_real_, qx.upper = NA_real_)
      df = rbind(df, aux); row.names(df) = NULL
    }

  }

  return(df[order(df$age), ])
}
