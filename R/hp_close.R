#' @title HP: Fitting the advanced ages of the life tables
#'
#' @description This function receives an object of the class `HP` fitted by the hp() function
#' and fits a closing method to expand the life tables dataset to a maximum age argument inputted
#' by the user.
#' There are four closing methods available: 'hp', 'plateau', 'linear', and 'gompertz'.
#' The 'linear' method can only be applied with HP objects following the lognormal variant of
#' the HP mortality law.
#'
#' @usage
#' hp_close (fit, method = c("hp", "plateau", "linear", "gompertz"),
#'  x0 = max(fit$data$x), max_age = 120, k = 7,
#'  weights = seq(from = 0, to = 1, length.out = 2*k+1),
#'  new_Ex = NULL, new_Dx = NULL)
#'
#' @param fit Object of the class `HP` fitted by the hp() function
#' @param method Character string specifying the closing method to be fitted, with them being: 'hp', 'plateau', 'linear' or 'gompertz'.
#' @param x0 Integer with the starting age the closing method will be fitted from. Default is the last age fitted by the 'HP' object.
#' @param max_age Integer with the maximum age the closing method will be fitted. Default age is '120'.
#' @param k Integer representing the size of the age interval to be mixed with the 'linear' or 'gompertz' closing methods for smooth graduation. If k = 0, no mixing will be applied.
#' @param weights Vector of weights of the closing method used in the mixture of the closing method and the fitted model made in the mixing age group. The vector's size should be equal to 2k+1. For a better understanding of this parameter and the mixture applied in this function, see Details.
#' @param new_Ex Vector with exposure of ages after the x0 input. This is an optional argument used in the 'linear' and 'gompertz' closing methods. If this argument is specified, then new_Dx also needs to be.
#' @param new_Dx Vector containing the death counts of the ages after the x0 input. This is also an optional argument used in the 'linear' and 'gompertz' closing methods. The length must be the same as new_Ex.
#'
#' @details
#' There are three types of age groups when the closing method is applied: a group
#' where only the HP-fitted model computes the death probabilities, followed by a
#' group in which the death probabilities are a mix (or more precise a weighted mean)
#' from the HP model and the closing method and followed by a group in which the
#' death probabilities are computed just by the closing method. The mix is applied
#' so the transition of the death probabilities of the ages between the fitted model
#' and the closing method occurs smoothly.
#'
#' The parameters 'x0' and 'k' define the mixing group age. The parameter 'x0'
#' indicates the center age of the group. The parameter 'k' is the range of ages
#' before 'x0' and after 'x0', so this group has a total of \eqn{2k + 1} age. Therefore,
#' the parameter 'weights' must have a length size equal to \eqn{2k + 1}. In this case,
#' the death probability is calculated as follows. Consider \eqn{model_x} and \eqn{close_x}
#' as the death probability of the fitted model and closing method in the age \eqn{x},
#' respectively. Then, the resulting death probability of the mix is calculated as:
#'
#' \eqn{q_x = w_x model_x + (1-w_x)close_x},
#'
#' where \eqn{w_x} represents the weight of the closing method in the age \eqn{x}.
#' This computation is applied to all elements in the MCMC chain of the fitted model,
#' resulting in a new chain of death probabilities. This procedure is applied only in
#' the linear and Gompertz methods.
#'
#' The four closing methods for life tables are:
#'
#' 1.'hp' method: Expands the previously adjusted HP model until the max_age argument.
#'
#' 2.'plateau' method: Keeps the death probability (qx) constant after the x0 argument.
#'
#' 3.'linear' method: Fits a linear regression starting at age x0 - k until the last age with data available (lognormal only).
#'
#' 4.'gompertz' method: Adopted as the closing method of the 2010-2012 English Life Table No. 17, fits the Gompertz mortality law via SIR using the same available data as the 'linear' method.
#'
#' @return Returns a `ClosedHP` class object with the predictive chains of the death probability
#' (qx) from first fitted age to max_age argument, the data utilized by the function and the
#' closing method chosen.
#'
#' @references Dodd, Erengul, Forster, Jonathan, Bijak, Jakub, & Smith, Peter 2018. “Smoothing mortality data: the English life table, 2010-12.” \emph{Journal of the Royal Statistical Society: Series A (Statistics in Society)}, 181(3), 717-735.
#'
#' @examples
#' ## Importing mortality data from the USA available on the Human Mortality Database (HMD):
#' data(USA)
#'
#' ## Selecting the exposure and the death count of the year 2010, ranging from 0 to 90 years old:
#' USA2010 = USA[USA$Year == 2010,]
#' x = 0:90
#' Ex = USA2010$Ex.Male[x+1]
#' Dx = USA2010$Dx.Male[x+1]
#'
#' ## Fitting a lognormal HP model:
#' fit = hp(x = x, Ex = Ex, Dx = Dx, model = "lognormal",
#'          M = 1000, bn = 0, thin = 10)
#'
#' ## Applying the closing function with different methods:
#' close1 = hp_close(fit, method = "hp", x0 = 90)
#' \donttest{close2 = hp_close(fit, method = "plateau", x0 = 90)
#' close3 = hp_close(fit, method = "linear", x0 = 80,
#'                   new_Ex = USA2010$Ex.Male[82:101],
#'                   new_Dx = USA2010$Dx.Male[82:101])
#' close4 = hp_close(fit, method = "gompertz", x0 = 70,
#'                   new_Ex = USA2010$Ex.Male[72:101],
#'                   new_Dx = USA2010$Dx.Male[72:101],
#'                   k = 5, max_age = 120)
#'
#' #### Using the other functions available in the package with the 'ClosedHP' object:
#'
#' ## qx estimation (See "?fitted.HP" in the BayesMortalityPlus package for more options):
#' fitted(close2)
#'
#' ## life expectancy (See "?expectancy.HP" for more options)
#' expectancy(close3, age = 0:110)
#'
#' ## plotting (See "?plot.HP" in the BayesMortalityPlus package for more options):
#' plot(close4)
#' g <- plot(list(close4, fit),
#'           colors = c("seagreen", "blue"),
#'           labels = c("Closed", "Model"))
#' # plotly::ggplotly(g)
#' }
#'
#' @seealso [fitted.HP()], [plot.HP()], [print.HP()] and [summary.HP()] for `ClosedHP` methods to native R functions [fitted()],
#'[plot()], [print()] and [summary()].
#'
#'[expectancy.HP()] and [Heatmap.HP()] for `ClosedHP` methods to compute and visualise the truncated life expectancy
#'via [expectancy()] and [Heatmap()] functions.
#'
#'
#' @include fun_aux.R
#' @include sir_gompertz.R
#'
#' @importFrom MASS mvrnorm
#'
#' @export
hp_close = function(fit, method = c("hp", "plateau", "linear", "gompertz"), x0 = max(fit$data$x),
                    max_age = 120, k = 7,  weights = seq(from = 0, to = 1, length.out = 2*k+1),
                    new_Ex = NULL, new_Dx = NULL){

  ## Pre-processing
  method = match.arg(method)

  ## Checklist
  if(!inherits(fit, "HP")) { stop("fit argument must be a 'HP' Object returned by hp() function.") }

  if (length(weights) != 2*k + 1) { stop("The length of the weights vector is not equal to 2k+1.") }

  if(x0 > max(fit$data$x)) { stop("x0 argument exceeds the maximum age of the model.") }

  if(x0 >= max_age) { stop("the choices of values for x0 and max_age are not consistent, x0 must be less than max_age.") }

  if(method %in% c("hp", "plateau")) { k = 0; weights = 0 }

  if(fit$info$model == "lognormal" & method %in% c("hp", "plateau")) { new_Ex = new_Dx = NULL }

  if (method == "linear" & fit$info$model != "lognormal"){ stop("linear closing method is only available for the lognormal model.") }

  if((method == "gompertz" | method == "linear") & ((!is.null(new_Ex) & is.null(new_Dx) | (is.null(new_Ex) & !is.null(new_Dx))))) { stop("gompertz and linear closing methods require new_Ex and new_Dx arguments.") }

  if(!is.null(new_Ex) & !is.null(new_Dx) & length(new_Ex) != length(new_Dx)) { stop("new_Ex and new_Dx lengths should be the same.") }


  ## Check if there are overlapping data between the model and the user input:
  min_age = min(fit$data$x)

  if(x0 < max(fit$data$x) & is.null(new_Ex)){
    new_Ex = fit$data$Ex[(x0+1-min_age):max(fit$data$x+1-min_age)]
    new_Dx = fit$data$Dx[(x0+1-min_age):max(fit$data$x+1-min_age)]
  }

  fit$data = fit$data[fit$data$x <= x0,]

  if(x0 - k < min_age) { stop("x0 or k arguments are not correct, they are not consistent with the initial age.") }

  ## Adding input data to the model data:
  age_last_data_Ex = x0 + length(new_Ex)
  age_last_data_Dx = x0 + length(new_Dx)
  new_Ex = c(fit$data$Ex, new_Ex)
  new_Dx = c(fit$data$Dx, new_Dx)

  if(max_age-age_last_data_Ex < 0) {max_age = age_last_data_Ex}

  ## Completing the data for the closing method:
  if(fit$info$model == "lognormal") {
    new_Ex = c(new_Ex, rep(NA_real_, max_age-age_last_data_Ex))
  }else{
    new_Ex = c(new_Ex, rep(new_Ex[length(new_Ex)], max_age-age_last_data_Ex))
  }
  new_Dx = c(new_Dx, rep(NA_real_, max_age-age_last_data_Dx))

  ## Data between 0 and the maximum age:
  full_Ex = c(new_Ex)
  full_Dx = c(new_Dx)

  ### CHECK
  if(method == "linear" | method == "gompertz"){
    data = data.frame(x = (x0-k):(age_last_data_Ex))
    data$Ex = full_Ex[data$x + 1 - min_age]
    data$Dx = full_Dx[data$x + 1 - min_age]
    data$qx = 1 - exp(-data$Dx/data$Ex)
    data$y = log(data$qx)

    if(nrow(data) < 2) { stop("Insufficient data to apply the closing method. Decrease the value of x0 argument or increase the value of k or try different data.") }

  }

  ## End length of the Markov chains:
  num_sim = nrow(fit$post.samples$mcmc_theta)

  ## Ages where the closing method will be applied:
  old_x = (x0 - k):max_age
  old_len = length(old_x)

  ## Matrix to save the fit:
  closed = matrix(0, nrow = num_sim, ncol = old_len)
  colnames(closed) = old_x

  ## Returns of the function: qx chains, x = min_age, ..., max_age
  ret = matrix(NA_real_, nrow = num_sim, ncol = max_age + 1 - min_age)

  ## Closing methods
  if (method == "hp"){

    # Fits the mortality curve to a group of parameters:
    if(fit$info$model == "lognormal"){
      for (i in 1:num_sim){
        hp = hp_curve(old_x, fit$post.samples$mcmc_theta[i, ])
        sim = exp(rnorm(old_len, log(hp), sqrt(fit$post.samples$sigma2[i])))
        closed[i, ] = sim/(1+sim)
      }
    }else if(fit$info$model == "binomial"){
      aux_Ex = full_Ex[old_x + 1]
      for (i in 1:num_sim){
        qx = 1 - exp(-hp_curve_9(old_x, fit$post.samples$mcmc_theta[i,]))
        qx = ifelse(qx > 1, 1, qx)
        qx = ifelse(qx < 0, 0, qx)
        sim = rbinom(old_len, trunc(aux_Ex), qx)
        closed[i, ] = sim/trunc(aux_Ex)
      }
    }else{
      aux_Ex = full_Ex[old_x + 1]
      for (i in 1:num_sim){
        qx = 1 - exp(-hp_curve_9(old_x, fit$post.samples$mcmc_theta[i,]))
        qx = ifelse(qx > 1, 1, qx)
        qx = ifelse(qx < 0, 0, qx)
        sim = rpois(old_len, aux_Ex*qx)
        closed[i, ] = sim/aux_Ex
      }
    }

  }else if(method == "plateau"){

    if(fit$info$model == "lognormal"){
      # for (i in 1:num_sim){
      #   hp = hp_curve(x0, fit$post.samples$mcmc_theta[i, ])
      #   sim = exp(rnorm(1, log(hp), sqrt(fit$post.samples$sigma2[i])))
      #   closed[i, ] = sim/(1+sim)
      # }
      # gets the death probability of x0 and applies it till max_age
      hp = hp_curve(x0, fit$post.samples$mcmc_theta)
      sim = exp(rnorm(num_sim, log(hp), sqrt(fit$post.samples$sigma2)))
      closed <- matrix(sim/(1+sim), num_sim, old_len)

    }else if(fit$info$model == "binomial"){
      aux_Ex = full_Ex[x0+1-min_age]
      # for (i in 1:num_sim){
      #   qx = 1 - exp(-hp_curve_9(x0, fit$post.samples$mcmc_theta[i,]))
      #   qx = ifelse(qx > 1, 1, qx)
      #   qx = ifelse(qx < 0, 0, qx)
      #   sim = rbinom(1, trunc(aux_Ex), qx)
      #   closed[i, ] = sim/trunc(aux_Ex)
      # }
      qx = 1 - exp(-hp_curve_9(x0, fit$post.samples$mcmc_theta))
      qx = ifelse(qx > 1, 1, qx)
      qx = ifelse(qx < 0, 0, qx)
      sim = rbinom(num_sim, rep(trunc(aux_Ex), num_sim), qx)
      closed <- matrix(sim/trunc(aux_Ex), num_sim, old_len)

    }else{
      aux_Ex = full_Ex[x0+1-min_age]
      # for (i in 1:num_sim){
      #   qx = 1 - exp(-hp_curve_9(x0, fit$post.samples$mcmc_theta[i,]))
      #   qx = ifelse(qx > 1, 1, qx)
      #   qx = ifelse(qx < 0, 0, qx)
      #   sim = rpois(1, aux_Ex*qx)
      #   closed[i, ] = sim/aux_Ex
      # }
      qx = 1 - exp(-hp_curve_9(x0, fit$post.samples$mcmc_theta))
      qx = ifelse(qx > 1, 1, qx)
      qx = ifelse(qx < 0, 0, qx)
      sim = rpois(num_sim, aux_Ex*qx)
      closed <- matrix(sim/aux_Ex, num_sim, old_len)
    }

  }else if(method == "linear"){

    mod = lm(y ~ x, data = data)
    pred = predict(mod, newdata = data.frame(x = old_x))

    X = model.matrix(mod)
    Xpred = cbind(1, old_x)
    C1 = t(X) %*% X
    Cpred = Xpred %*% chol2inv(chol(C1)) %*% t(Xpred)
    RMAT = (diag(old_len) + Cpred)

    for (i in 1:num_sim){
      sig = sqrt(fit$post.samples$sigma2[i])
      SIGMApred = sig * RMAT

      sim_vals = MASS::mvrnorm(1, mu = pred, Sigma = SIGMApred)
      closed[i, ] = exp(sim_vals)
    }

  }else if(method == "gompertz"){

    param = sir_gompertz(fit, data, resampling_size = num_sim)

    if(fit$info$model == "lognormal"){

      for (i in 1:num_sim){
        gomp = param[i,1]*exp(param[i,2]*old_x)
        sim = exp(rnorm(old_len, log(gomp), sqrt(fit$post.samples$sigma2[i])))
        closed[i, ] = sim/(1+sim)
      }
    }else if(fit$info$model == "binomial"){
      aux_Ex = full_Ex[old_x + 1]
      for (i in 1:num_sim){
        qx = 1 - exp(-param[i,1]*exp(param[i,2]*old_x))
        qx = ifelse(qx > 1, 1, qx)
        qx = ifelse(qx < 0, 0, qx)
        sim = rbinom(old_len, trunc(aux_Ex), qx)
        closed[i, ] = sim/trunc(aux_Ex)
      }
    }else{
      aux_Ex = full_Ex[old_x + 1]
      for (i in 1:num_sim){
        qx = 1 - exp(-param[i,1]*exp(param[i,2]*old_x))
        qx = ifelse(qx > 1, 1, qx)
        qx = ifelse(qx < 0, 0, qx)
        sim = rpois(old_len, aux_Ex*qx)
        closed[i, ] = sim/aux_Ex
      }
    }
  }

  ## qx error margin (default was 0.02)
  eps = 0.01
  # Preventing death probabilities above 1:
  closed = apply(closed, 2, function(x) ifelse(x < 1 - eps, x, 1))

  new_age = min_age:max_age
  new_len = length(new_age)
  fitted = matrix(0, nrow = num_sim, ncol = new_len)
  colnames(fitted) = new_age

  ## Checking the model used:
  if(fit$info$model == "lognormal"){
    for (i in 1:num_sim){
      hp = hp_curve(new_age, fit$post.samples$mcmc_theta[i, ])
      sim = exp(rnorm(new_len, log(hp), sqrt(fit$post.samples$sigma2[i])))
      fitted[i, ] = sim/(1+sim)
    }
  }else if(fit$info$model == "binomial"){
    for (i in 1:num_sim){
      qx = 1 - exp(-hp_curve_9(new_age, fit$post.samples$mcmc_theta[i,]))
      qx = ifelse(qx > 1, 1, qx)
      qx = ifelse(qx < 0, 0, qx)
      sim = rbinom(new_len, trunc(full_Ex), qx)
      fitted[i, ] = sim/trunc(full_Ex)
    }
  }else{
    for (i in 1:num_sim){
      qx = 1 - exp(-hp_curve_9(new_age, fit$post.samples$mcmc_theta[i,]))
      qx = ifelse(qx > 1, 1, qx)
      qx = ifelse(qx < 0, 0, qx)
      sim = rpois(new_len, full_Ex*qx)
      fitted[i, ] = sim/full_Ex
    }
  }

  # Model only indexes: age 0 till x0 - k - 1
  idx_mod_only = min_age:(x0 - k - 1) + 1 - min_age
  # Mix indexes: age x0 - k till x0 + k
  idx_mix = (x0 - k):(x0 + k) + 1 - min_age
  # Closing method only indexes: age x0 + k + 1 till max_age
  idx_close = (x0 + k + 1):max_age + 1 - min_age

  idx_mnc = c(idx_mix, idx_close)

  ret[ , idx_mod_only] = fitted[ , idx_mod_only]
  ret[ , idx_mnc] = closed

  # Mix
  for (i in 1:num_sim) {
    ret[i, idx_mix] = weights * ret[i, idx_mix] + (1 - weights) * fitted[i, idx_mix]
  }

  return(structure(list(qx = ret,
                        data = data.frame(x = new_age,
                                          Ex = full_Ex,
                                          Dx = full_Dx,
                                          qx = 1-exp(-full_Dx/full_Ex)),
                        method = method), class = "ClosedHP"))
}
