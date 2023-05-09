#' @title Predictive intervals for deaths probabilities.
#'
#' @description
#' This function returns the predictive intervals for death probabilities of a
#' mortality graduation provide by hp(), dlm(), hp_close() or dlm_close() functions
#'
#' @usage
#' qx_ci(fit, age = NULL, Ex = NULL, prob = 0.95)
#'
#' @param fit Object of the following classes: 'HP', 'DLM', 'ClosedHP', or 'ClosedDLM'.
#' @param age Vector with the ages to compute the predictive intervals.
#' @param Ex Vector with the exposures of the selected ages. Its length must be equal to the age vector. This argument is only necessary when using the Poisson and the Binomial distributions.
#' @param prob Coverage probability of the predictive intervals.
#'
#' @return A data.frame object with the selected ages and the respective predictive intervals of the death probabilities.
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
#' ## Fitting poisson HP model and DLM:
#' fit = hp(x = x, Ex = Ex, Dx = Dx, model = "poisson",
#'          M = 1000, bn = 0, thin = 10)
#'
#' ## predictive intervals
#' qx_ci(fit)
#' qx_ci(fit, age = 0:110, Ex = USA2000$Ex.Total, prob=0.99)
#'
#' @include fun_aux.R
#'
#' @export
qx_ci <- function(fit, age = NULL, Ex = NULL, prob = 0.95){

  if(inherits(fit, "HP")){
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
      fitted = matrix(NA, nrow = nrow(fit$post.samples$mcmc_theta), ncol = length(age))

      for (i in 1:nrow(fitted)){
        qx = 1 - exp(-hp_curve_9(age, fit$post.samples$mcmc_theta[i,]))
        qx = ifelse((qx < 0 | qx > 1), NA, qx)
        sim = rbinom(length(age), trunc(Ex), qx)
        fitted[i,] = sim/Ex
      }
    }else if(fit$info$model == "poisson"){

      age_out = age[is.na(Ex)]
      age = age[!is.na(Ex)] ## Removing ages with no exposures
      Ex = Ex[!is.na(Ex)] ## Removing exposures with NA values
      fitted = matrix(NA, nrow = nrow(fit$post.samples$mcmc_theta), ncol = length(age))

      for(i in 1:nrow(fitted)){
        qx = 1 - exp(-hp_curve_9(age, fit$post.samples$mcmc_theta[i,]))
        qx = ifelse((qx < 0 | qx > 1), NA, qx)
        sim = rpois(length(age), lambda = Ex*qx)
        fitted[i,] = sim/Ex
      }
    }else{

      age_out = NULL
      fitted = matrix(NA, nrow = nrow(fit$post.samples$mcmc_theta), ncol = length(age))

      for(i in 1:nrow(fitted)){
        hp <- hp_curve(age, fit$post.samples$mcmc_theta[i,])
        sim = rnorm(length(age), log(hp), sqrt(fit$post.samples$sigma2[i]))
        fitted[i,] <- exp(sim)/(1+exp(sim))
      }
    }


    qi = apply(fitted, 2, quantile, (1-prob)/2, na.rm = T)
    qs = apply(fitted, 2, quantile, (1+prob)/2, na.rm = T)

    aux = data.frame(age = age, qi = qi, qs = qs)
    aux[!(aux$qi > 0), 2] = 0
    aux[!(aux$qs < 1), 3] = 1

    if(length(age_out) > 0){
      aux2 <- data.frame(age = age_out, qi = NA, qs = NA)
      aux <- rbind(aux, aux2)
      aux <- aux[order(aux$age),]
    }

    return(aux)

  }else if(inherits(fit, "ClosedHP")){
    if(fit$method == "Mix"){
      return(data.frame(age = fit$data$x, qi = NA, qs = NA))
    }

    fitted = fit$qx
    close_age = fit$data$x

    qi = apply(fitted, 2, quantile, (1-prob)/2, na.rm = T)
    qs = apply(fitted, 2, quantile, (1+prob)/2, na.rm = T)

    df = data.frame(age = close_age, qi = qi, qs = qs)
    df[!(df$qi > 0), 2] = 0
    df[!(df$qs < 1), 3] = 1

    if(!is.null(age)){

      df = df[(close_age %in% age), ]

      if(any(!(age %in% close_age))){
        age_not_fitted = age[!(age %in% close_age)]
        aux = data.frame(age = age_not_fitted, qi = NA, qs = NA)
        df = rbind(df, aux); row.names(df) = NULL
      }

    }

    return(df[order(df$age), ])

  }else if(inherits(fit,"DLM")){

    t = ncol(fit$mu)
    n = nrow(fit$mu)

    fitted = matrix(NA, nrow = n, ncol = t)

    for(i in 1:n){
      sim = rnorm(t, fit$mu[i,], sqrt(fit$sig2[i]))
      fitted[i,] <- exp(sim)
    }

    qi = apply(fitted, 2, quantile, (1-prob)/2, na.rm = T)
    qs = apply(fitted, 2, quantile, (1+prob)/2, na.rm = T)

    aux = data.frame(age = fit$info$ages, qi = qi, qs = qs)
    aux[!(aux$qi > 0), 2] = 0
    aux[!(aux$qs < 1), 3] = 1

    if(!is.null(age)) aux = aux[(aux$age %in% age), ]

    return(aux)

  }else if(inherits(fit,"ClosedDLM")){

    fitted = fit$qx
    close_age = fit$info$ages

    qi = apply(fitted, 2, quantile, (1-prob)/2, na.rm = T)
    qs = apply(fitted, 2, quantile, (1+prob)/2, na.rm = T)

    df = data.frame(age = close_age, qi = qi, qs = qs)
    # df[!(df$qi > 0), 2] = 0
    # df[!(df$qs < 1), 3] = 1

    if(!is.null(age)) df = df[(df$age %in% age), ]

    return(df[order(df$age), ])

  }
}
