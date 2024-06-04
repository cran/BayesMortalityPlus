#' @title DLM: Fitting the advanced ages of the life tables
#'
#' @description This function receives an object of the class `DLM` fitted by the dlm() function
#' and fits a closing method to expand the life tables dataset to a maximum age argument inputed
#' by the user.
#' There are three closing methods available: 'linear', 'gompertz' and 'plateau'.
#'
#' @usage
#' dlm_close(fit, method = c("linear", "gompertz", "plateau"),
#'           x0 = max(fit$info$ages), max_age = 120, k = 7,
#'           weights = seq(from = 0, to = 1, length.out = k),
#'           new_data = NULL)
#'
#' @param fit Object of the class `DLM` adjusted by the dlm() function.
#' @param method Character string specifying the closing method to be fitted, with them being: 'plateau', 'linear' or 'gompertz'.
#' @param x0 Integer with the starting age the closing method will be fitted from. Default is the last age fitted by the 'DLM' object.
#' @param max_age Integer with the maximum age the closing method will be fitted. Default age is '120'.
#' @param k Integer representing the size of the age-interval to be mixed with the 'linear' or 'gompertz' closing methods for a smooth graduation. If k = 0, no mixing will be made. Default: 7.
#' @param weights Vector of weights of the closing method used in the mixture of the closing method and the fitted model made in the mixing age group. The vector's size should be equal to 2k+1. For a better understanding of this parameter and the mixture applied in this function, see Details.
#' @param new_data Vector containing the log mortality rates of ages after the x0 input. This is an optional argument used in the 'linear' and 'Gompertz' closing methods.
#'
#' @details
#' #' There are three types of age groups when the closing method is applied: a group
#' where only the fitted model (DLM) computes the death probabilities, followed by a
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
#' where \eqn{w_x} represents the weight of the closing method in the age \eqn{x}. This
#' procedure is applied only in the linear and Gompertz methods.
#'
#'
#' The three closing methods implemented by the function are:
#'
#' 1.'linear' method: Fits a linear regression starting at age x0 - k until the last age with data available
#'
#' 2.'gompertz' method: Used as the closing method of the 2010-2012 English Life Table No. 17, fits the Gompertz mortality law via SIR using the same available data as the 'linear' method.
#'
#' 3.'plateau' method: Keeps the death probability (qx) constant after the x0 argument.
#'
#' @return Returns a `ClosedDLM` class object with the predictive chains of the death probability
#' (qx) from first fitted age to max_age argument, the data information utilized by the function and the
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
#' x = 0:100
#' Ex = USA2010$Ex.Male[x+1]
#' Dx = USA2010$Dx.Male[x+1]
#' y <- log(Dx/Ex)
#'
#' fit <- dlm(y, M = 100)
#'
#' ## Applying the closing function with different methods:
#' close1 = dlm_close(fit, method = "plateau")
#'
#' ### Getting new data for the linear and gompertz methods:::
#' x2 = 101:110
#' Ex2 = USA2010$Ex.Male[x2+1]
#' Dx2 = USA2010$Dx.Male[x2+1]
#' y2 <- log(Dx2/Ex2)
#'
#' close2 = dlm_close(fit, method = "linear",
#'                   new_data = y2)
#'
#' #### Using the other functions available in the package with the 'ClosedDLM' object:
#'
#' ## qx estimation (See "?fitted" in the BayesMortalityPlus package for more options):
#' fitted(close2)
#'
#' ## life expectancy (See "?expectancy.DLM" for more options)
#' expectancy(close2, age = seq(0,120,by=20), graph = FALSE)
#'
#' ## plotting (See "?plot" in the BayesMortalityPlus package for more options):
#' plot(list(close1, close2, fit),
#'      colors = c("red4","seagreen", "blue"),
#'      labels = c("Plateau method","Linear method", "DLM fitted"),
#'      plotData = FALSE)
#'
#'
#' @include gompertz_dlm.R
#'
#' @importFrom MASS mvrnorm
#'
#' @seealso [fitted.DLM()], [plot.DLM()], [print.DLM()] and [summary.DLM()] for `ClosedDLM` methods to native R functions [fitted()],
#'[plot()], [print()] and [summary()].
#'
#'[expectancy.DLM()] and [Heatmap.DLM()] for `ClosedDLM` methods to compute and visualise the truncated life expectancy
#'via [expectancy()] and [Heatmap()] functions.
#'
#'
#' @export
dlm_close = function(fit, method = c("linear", "gompertz", "plateau"), x0 = max(fit$info$ages),
                     max_age = 120, k = 7,  weights = seq(from = 0, to = 1, length.out = k),
                     new_data = NULL){

  ## Pre-processing
  method = match.arg(method)

  ## Checklist
  if(!inherits(fit, "DLM")) { stop("fit argument must be a 'DLM' Object returned by dlm() function.") }

  if (length(weights) != k) { stop("The length of the weights vector is not equal to k.") }

  if(x0 > max(fit$info$ages)) { stop("x0 argument exceeds the maximum age of the model.") }

  if(x0 >= max_age) { stop("The choices of values for x0 and max_age are not consistent, x0 must be less than max_age.") }

  if(method == "plateau") { k = 1; weights = 0; new_data = NULL }

  #if((method == "gompertz" | method == "linear") & is.null(new_data)) { stop("gompertz and linear closing methods require new_data argument.") }


  ## Check if there are overlapping data between the model and the user input:
  min_age = min(fit$info$ages)

  if((x0 < max(fit$info$ages)) & is.null(new_data)){
    new_data = fit$info$y[(x0+1-min_age):max(fit$info$ages+1-min_age)]
  }

  fit$info$y = fit$info$y[fit$info$ages <= x0]

  if(x0 - k < min_age) { stop("x0 or k arguments are not correct, they are not consistent with the initial age.") }

  ## Adding input data to the model data:
  age_last_data = x0 + length(new_data)
  new_data = c(fit$info$y, new_data)

  if(max_age-age_last_data < 0) {max_age = age_last_data}

  ## Completing the data for the closing method:
  new_data = c(new_data, rep(NA_real_, max_age-age_last_data))

  ## Checking if the model starts at the age 0:
  if(min_age != 0){ pre_data = rep(NA_real_, length(1:min_age)) }else{ pre_data = NULL }

  ## Data between 0 and the maximum age:
  full_data = c(pre_data, new_data)

  if(method == "linear" | method == "gompertz"){
    data = data.frame(x = (x0-k+1):(age_last_data))
    data$y = full_data[data$x + 1]

    if(nrow(data) < 2) { stop("Insufficient data to apply the closing method. Decrease the value of x0 argument or increase the value of k or try different data.") }
  }

  ## End length of the Markov chains:
  num_sim = length(fit$sig2)

  ## Ages where the closing method will be applied:
  old_x = (x0 - k + 1):max_age
  old_len = length(old_x)

  ## Matrix to save the fit:
  closed = matrix(0, nrow = num_sim, ncol = old_len)
  colnames(closed) = old_x

  ## Returns of the function: qx chains, x = min_age, ..., max_age
  ret = matrix(NA_real_, nrow = num_sim, ncol = max_age + 1)

  ## Closing methods
  if(method == "plateau"){

    #gets the death probability of x0 and applies it till max_age
    sim <- 1 - exp(-exp(rnorm(num_sim, fit$mu[,x0-min_age+1], sqrt(0.01*fit$sig2))))
    closed <- matrix(sim, num_sim, old_len)
    colnames(closed) = old_x

    # for (i in 1:num_sim){
    #   sim = rnorm(1, fit$mu[i,x0-min_age+1], sqrt(fit$sig2[i]))
    #   closed[i, ] = exp(sim)
    # }

  }else if(method == "linear"){

    mod = lm(y ~ x, data = data)
    pred = predict(mod, newdata = data.frame(x = old_x))

    X = model.matrix(mod)
    Xpred = cbind(1, old_x)
    C1 = t(X) %*% X
    Cpred = Xpred %*% chol2inv(chol(C1)) %*% t(Xpred)
    RMAT = (diag(old_len) + Cpred)

    for (i in 1:num_sim){
      sig = sqrt(0.01*fit$sig2[i])
      SIGMApred = sig * RMAT

      sim_vals = MASS::mvrnorm(1, mu = pred, Sigma = SIGMApred)
      closed[i, ] = 1 - exp(-exp(sim_vals))
    }

  ####################################################################################################
  }else if(method == "gompertz"){

    param = sir_gompertz_dlm(fit, data, resampling_size = num_sim)

    for (i in 1:num_sim){
      gomp = param[i,1]*exp(param[i,2]*old_x)
      sim = 1 - exp(-exp(rnorm(old_len, log(gomp), sqrt(0.01*fit$sig2[i]))))
      closed[i, ] = sim
    }
  }
  ####################################################################################################
  ## qx error margin (default was 0.02)
  eps = 0.01
  # Preventing death probabilities above 1:
  closed = apply(closed, 2, function(x) ifelse(x < 1 - eps, x, 1))

  new_age = min(fit$info$ages):max_age
  new_len = length(new_age)
  fitted = matrix(NA_real_, nrow = num_sim, ncol = new_len)
  colnames(fitted) = new_age

  ####################################################################################################
  t = ncol(fit$mu)
  for (i in 1:num_sim){
    sim = rnorm(t, fit$mu[i,], sqrt(0.01*fit$sig2[i]))
    fitted[i, (min_age+1):(t+min_age)] = 1 - exp(-exp(sim))
  }
  ####################################################################################################
  # Model only indexes: age 0 till x0 - k - 1
  idx_mod_only = 0:(x0 - k) + 1
  # Mix indexes: age x0 - k till x0 + k
  idx_mix = (x0 - k + 1):x0 + 1
  # Closing method only indexes: age x0 + k + 1 till max_age
  idx_close = (x0 + 1):max_age + 1

  idx_mnc = c(idx_mix, idx_close)

  ret[ , idx_mod_only] = fitted[ , idx_mod_only]
  ret[ , idx_mnc] = closed

  # Mix
  for (i in 1:num_sim) {
    ret[i, idx_mix] = weights * ret[i, idx_mix] + (1 - weights) * fitted[i, idx_mix]
  }

  return(structure(list(qx = ret[, (min(fit$info$ages):max_age + 1)],
                        info = list(ages = new_age,
                                    y = full_data),
                        method = method), class = "ClosedDLM"))
}
