#' @title HP Model mixture
#'
#' @description This function mixes the fitted mortality table of the HP model with another mortality
#' table provided by the user.
#'
#' @usage
#' hp_mix (fit, mu_post, weights = NULL, mix_age,
#'  x0_prior, x0_post, max_age)
#'
#' @param fit Object of the class 'HP' fitted by the hp() function.
#' @param mu_post Vector with mortality rates considered in the mix.
#' @param weights Positive vector specifying the weights considered in the mix.
#' @param mix_age Positive vector specifying the age range in the mixture.
#' @param x0_prior Non-negative number indicating the initial age of the fitted HP model.
#' @param x0_post Non-negative number indicating the initial age of the mortality table provided by the user.
#' @param max_age Positive number indicating the final age in the mixture.
#'
#' @return Return the posterior distribution for qx.
#'
#' @examples
#' ## Importing mortality data from the USA available on the Human Mortality Database (HMD):
#' data(USA)
#'
#' ## Selecting the exposure and death count of the 2010 and 2013 male populations ranging
#' ## from 0 to 90 years old
#' USA2010 = USA[USA$Year == 2010,]
#' x = 0:90
#' Ex = USA2010$Ex.Male[x+1]
#' Dx = USA2010$Dx.Male[x+1]
#'
#' USA2013 = USA[USA$Year == 2013,]
#' Ex2 = USA2013$Ex.Male[x+1]
#' Dx2 = USA2013$Dx.Male[x+1]
#'
#' ## Fitting HP model for 2010 data and calculating the mortality rates of 2013
#' fit = hp(x = x, Ex = Ex, Dx = Dx,
#'          M = 1000, bn = 0, thin = 10)
#' tx_2013 = 1 - exp(-Dx2/Ex2)
#'
#' ## Mixing fitted model and mortality rates of 2013:
#' mix <- hp_mix(fit, tx_2013, x0_prior = 0, x0_post = 0, mix_age = c(50,90),
#'               max_age = 90)
#'
#' ## Obtaining the new estimated mortality table (after mixture):
#' qx_mix<- apply(mix$qx, 2, median, na.rm = TRUE)
#' qx_mix
#'
#' @include fun_aux.R
#'
#' @export
hp_mix <- function(fit, mu_post, weights = NULL, mix_age, x0_prior, x0_post, max_age) {

  if(!inherits(fit, "HP")) stop("fit must be an object of the class HP.")
  if(x0_prior < 0) stop("x0_prior is the initial age used in the fitted HP model.")
  if(x0_post < 0) stop("x0_post is the initial age of the mortality rates used in the mixture.")
  if(mix_age[1] > mix_age[2]) { aux <- mix_age; mix_age[1] <- aux[2]; mix_age[2] <- aux[1]; rm(aux) }
  if(max_age < x0_prior) stop("max_age must be bigger than x0_prior.")
  if(max_age < x0_post) stop("max_age must be bigger than x0_post.")

  ## Verificar se weights foi passado corretamente
  if(is.null(weights)){
    prior_weights = seq(from = 1, to = 0, length.out = mix_age[2] - mix_age[1] + 1)
  }else{
    if(any(weights < 0)){
      warning("weights cannot have negative values.")
      prior_weights = seq(from = 1, to = 0, length.out = mix_age[2] - mix_age[1] + 1)
    }else if(any(weights > 1)){
      if(length(weights) != length(mix_age[2] - mix_age[1] + 1)) stop("Length of the vector 'weights' must be equal to the range of mix_age.")
      prior_weights = weights/max(weights)
    }else{
      if(length(weights) != (mix_age[2] - mix_age[1] + 1)) stop("Length of the vector 'weights' must be equal to the range of mix_age.")
      prior_weights = weights
    }
  }

  posterior_weights = 1 - prior_weights

  ## calcula qx das cadeias do fit
  age_fitted <- x0_prior:max_age
  mu_prior = matrix(NA, nrow = nrow(fit$post.samples$mcmc_theta), ncol = length(age_fitted))
  for (i in 1:nrow(mu_prior)){
    mu_prior[i,] = 1 - exp(-hp_curve_9(age_fitted, fit$post.samples$mcmc_theta[i,]))
    mu_prior[i,] = ifelse((mu_prior[i,] < 0 | mu_prior[i,] > 1), NA, mu_prior[i,])
  }

  ## Mistura
  mix_int = mix_age[1]:mix_age[2] + 1
  for (i in 1:nrow(mu_prior)) {
    mu_prior[i,mix_int-x0_prior] = prior_weights*mu_prior[i,mix_int-x0_prior] + posterior_weights*mu_post[mix_int-x0_post]
  }


  return(structure(list(qx = mu_prior,
                        data = data.frame(x = age_fitted,
                                          Ex = fit$data$Ex[age_fitted - min(fit$data$x)+1],
                                          Dx = fit$data$Dx[age_fitted - min(fit$data$x)+1],
                                          qx = 1-exp(-fit$data$Dx[age_fitted - min(fit$data$x)+1]/fit$data$Ex[age_fitted - min(fit$data$x)+1])),
                        info = list(mix_age = mix_age, x0_prior = x0_prior, x0_post = x0_post, max_age = max_age, mu_post = mu_post),
                        method = "Mix"),
                   class = "ClosedHP"))
}
