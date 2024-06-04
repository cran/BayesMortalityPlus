#' @title Dynamic Linear Model for mortality table graduation
#'
#' @description
#' This function fits a Dynamic Linear Model (DLM) for mortality data following
#' a Bayesian framework using Forward Filtering Backward Sampling algorithm to compute the posterior distribution.
#' The response variable is the log of the mortality rate, and it is modeled specifying the matrices Ft and Gt from the DLM equations.
#' Furthermore, the discount factor is used to control the smoothness of the fitted model. By default, a
#' linear growth model is specified.
#'
#' @usage
#' dlm(y, Ft = matrix(c(1,0), nrow = 1), Gt = matrix(c(1,0,1,1), 2), delta = 0.85,
#'     prior = list(m0 = rep(0, nrow(Gt)), C0 = diag(100, nrow(Gt))),
#'     prior.sig2 = list(a = 0, b = 0), M = 2000, ages = 0:(length(y)-1))
#'
#' @param y Numeric vector of log mortality rates.
#' @param Ft 1xp Matrix that specifies the observation equation, where p is the number of parameters. By default, 'Ft = matrix(c(1,0), nrow = 1)'.
#' @param Gt pxp Matrix that specifies the system equations. By default, Gt = matrix(c(1,0,1,1), 2).
#' @param delta Positive real value or real vector of the same length as y with values in the '(0, 1)' interval specifying the discount factor for each age. A higher value of delta results in a higher smoothness of the fitted curve. If a single value is defined, this same value is used for all ages. By default, delta is '0.85'.
#' @param prior A list with the prior mean vector \eqn{(m_0)} and covariance matrix \eqn{(C_0)} of \eqn{\theta_0} (state vector at time (age) t = 0). By default mean of zeros and diagonal matrix with a common variance 100 is used. Each element of the list must be named accordingly with the parameter (m0 for mean vector and C0 for covariance matrix).
#' @param prior.sig2 A list with the prior parameters (a, b) of Inverted Gamma distribution for \eqn{\sigma^2}. Each element of the list must be named accordingly with the parameter (a for shape parameter and b for scale parameter).
#' @param M Positive integer that indicates the sampling size from the posterior distributions. The default value is 2000.
#' @param ages Numeric vector of the ages fitted. Default is '0:(length(y)-1)'.
#'
#' @details
#' Let \eqn{Y_t} be the log mortality rate at age \eqn{t}. A DLM is specified as follows:
#'
#' For \eqn{t = 0}:
#'
#'  \eqn{\theta_0 \sim N_p (m_0, C_0)}
#'
#' Now, for \eqn{t \geq 1}:
#'
#' The observation equation:
#'
#'  \eqn{Y_t = F_t \theta_t + v_t}
#'
#' The system equation:
#'
#'  \eqn{\theta_t = G_t \theta_{t-1} + w_t}
#'
#' Where \eqn{F_t} and \eqn{G_t} are known matrices. \eqn{v_t} and \eqn{w_t} are independent
#' random errors with \eqn{v_t \sim N(0, \sigma^2)} and \eqn{w_t \sim N(0, \sigma^2 W_t)}. We
#' use the discount factors \eqn{\delta} to specify \eqn{W_t} as \eqn{W_t = C_t(1-\delta)/\delta},
#' where \eqn{C_t} is the conditional covariance matrix of \eqn{\theta_t}. So, if
#' \eqn{\delta = 0} there is no loss information as \eqn{t} increase (completely reducing the
#' smoothness of the fitted curve). \eqn{\delta} can be specified as a single value for all ages
#' or as a vector in which each element is associated with an age.
#'
#' A scheme described by (Petris et al, 2009) for conjugated inference is used.
#' For more details, see (Petris et al, 2009).
#'
#' @return A DLM class object.
#' \item{mu}{Posterior samples from \eqn{\mu_t = F_t \theta_t}, for all t.}
#' \item{theta}{Posterior samples from \eqn{\theta_t}, for all t.}
#' \item{sig2}{Posterior samples from \eqn{\sigma^2}.}
#' \item{param}{A list with the states parameters for filtering distribution (mt, Ct), predictive distribution (ft, Qt), smoothing distribution (as, Rs), and parameters of the posterior distribution for variance (alpha, beta).}
#' \item{info}{A list with some informations of the fitted model: the specification of \eqn{F_t} and \eqn{G_t} matrices, the data y and the ages, the discount factor \eqn{delta} value specified and priors informations.}
#'
#' @references Campagnoli, P., Petris, G., and Petrone, S. (2009). \emph{Dynamic linear models with R}. Springer-Verlag New York.
#'
#' @examples
#' ## Importing mortality data from the USA available on the Human Mortality Database (HMD):
#' data(USA)
#'
#' ## Selecting the log mortality rate of the 2010 male population ranging from 0 to 100 years old
#' USA2010 = USA[USA$Year == 2010,]
#' x = 0:100
#' Ex = USA2010$Ex.Male[x+1]
#' Dx = USA2010$Dx.Male[x+1]
#' y = log(Dx/Ex)
#'
#' ## Fitting DLM
#' fit = dlm(y)
#' print(fit)
#' summary(fit)
#'
#' ## Using other functions available in the package:
#' ## plotting (See "?plot.DLM" in the BayesMortality package for more options):
#' plot(fit)
#'
#' ## qx estimation (See "?fitted.DLM" in the BayesMortality package for more options):
#' fitted(fit)
#'
#' ## chain's plot (See "?plot_chain" for more options):
#' plot_chain(fit, param = c("mu[0]", "mu[100]"))
#'
#' ## Varying discount factor
#' fit2 = dlm(y, delta = c(rep(0.8, 36), rep(0.9, 65)))
#' plot(fit2)
#'
#'
#' @include ffbs.R
#'
#' @importFrom mvtnorm rmvt
#'
#'@seealso [fitted.DLM()], [predict.DLM()], [plot.DLM()], [print.DLM()] and [summary.DLM()] for `DLM` methods to native R functions [fitted()],
#'[plot()], [print()] and [summary()].
#'
#'[expectancy.DLM()] and [Heatmap.DLM()] for `DLM` methods to compute and visualise the truncated life expectancy
#'via [expectancy()] and [Heatmap()] functions.
#'
#'[dlm_close()] for close methods to expand the life tables.
#'
#'[plot_chain()] to visualise the markov chains, respectively.
#'
#' @export
dlm <- function(y, Ft = matrix(c(1,0), nrow = 1), Gt = matrix(c(1,0,1,1), 2), delta = 0.85,
                prior = list(m0 = rep(0, nrow(Gt)), C0 = diag(100, nrow(Gt))),
                prior.sig2 = list(a = 0, b = 0), M = 2000, ages = 0:(length(y)-1)){

  ## Validation
  if(is.vector(Ft)) {Ft = t(as.matrix(Ft))}
  if(nrow(Ft) != 1) stop("Ft must be a matrix with the following dimensions: 1 row and p columns.")
  if(!(is.matrix(Gt))) {Gt = as.matrix(Gt)}
  if(ncol(Ft) != nrow(Gt)) stop("Matrices Ft and Gt are not well defined.")
  if(ncol(Gt) != nrow(Gt)) stop("Gt must be a square matrix.")
  if(length(prior$m0) != nrow(Gt)) stop("Dimension of prior mean does not match the dimension of matrix Gt.")
  if(nrow(prior$C) != nrow(Gt)) stop("Dimension of prior covariance matrix does not match the dimension of matrix Gt.")
  if(ncol(prior$C) != nrow(Gt)) stop("Dimension of prior covariance matrix does not match the dimension of matrix Gt.")
  if(!(length(delta) == length(y)) & length(delta) != 1) stop("delta must be the same length of y")
  if(any(delta <= 0 | delta > 1)) stop("delta must be in the interval (0,1].")

  ## Auxiliary objects
  t = length(y)
  p = length(prior$m0)
  fit <- list()

  ## Filtering
  filter = ff(y = y, Ft = Ft, Gt = Gt, m0 = prior$m0, C0 = prior$C0, delta = delta,
              alpha0 = prior.sig2$a, beta0 = prior.sig2$b)

  ## Smoothing
  smooth = bs(m = filter$m, C = filter$C, a = filter$a, R = filter$R, Gt = Gt,
              alpha = filter$alpha, beta = filter$beta)

  ## Sampling
  sig2 = 1/rgamma(M, smooth$alpha, smooth$beta)
  theta <- array(NA, dim = c(M, t, p))
  mu <- matrix(NA, nrow = M, ncol = t)

  if(p == 1){
    for(i in 1:t){
      theta[,i,] = aux = rt(M, df = 2*smooth$alpha)*sqrt(c(smooth$Rs[i,,])*(smooth$beta/smooth$alpha)) + smooth$as[i,]
      mu[,i] <- c(aux%*%t(Ft))
    }
  }else{
    for(i in 1:t){
      theta[,i,] = aux = mvtnorm::rmvt(M, sigma = smooth$Rs[i,,]*(smooth$beta/smooth$alpha),
                                       delta = smooth$as[i,], df = 2*smooth$alpha, type = "shifted")
      mu[,i] <- c(aux%*%t(Ft))
    }
  }


  fit$mu = mu
  fit$theta = theta
  fit$sig2 = sig2
  # fit$Wt = filter$Wt
  fit$param = list(mt = filter$m, Ct = filter$C,
                   ft = filter$f, Qt = filter$Q,
                   as = smooth$as, Rs = filter$Rs,
                   alpha = smooth$alpha, beta = smooth$beta)
  fit$info = list(y = y,
                  ages = ages,
                  Ft = Ft,
                  Gt = Gt,
                  delta = delta,
                  prior = prior,
                  prior.sig2 = prior.sig2)

  return(structure(fit, class = "DLM"))
}
