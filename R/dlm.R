#' @title Dynamic Linear Model for mortality table graduation
#'
#' @description
#' This function fits a Dynamic Linear Model (DLM) for mortality data following
#' a Bayesian framework using Forward Filtering Backward Sampling algorithm and Markov chain Monte Carlo
#' Gibbs sampler to compute the posterior distribution. The response variable is the log of mortality rate,
#' and it is modeled specifying the matrices Ft and Gt from the DLM equations. Furthermore,
#' the discount factor is used to control the smoothness of the fitted model. By default, a
#' linear growth model is specified.
#'
#' @usage
#' dlm(y, Ft = matrix(c(1,0), nrow = 1), Gt = matrix(c(1,0,1,1), 2),
#'  delta = 0.85, prior = list(m0 = rep(0, nrow(Gt)), C0 = diag(100, nrow(Gt))),
#'  prior.sig2 = list(a = 0.01, b = 0.01), M = 5000, bn = 3000, thin = 1,
#'  ages = 0:(length(y)-1))
#'
#' @param y Numeric vector of log mortality rates.
#' @param Ft 1xp Matrix that specifies the observation equation, where p is the number of parameters. By default, 'Ft = matrix(c(1,0), nrow = 1)'.
#' @param Gt pxp Matrix that specifies the system equations. By default, Gt = matrix(c(1,0,1,1), 2).
#' @param delta Positive number in '(0, 1)' interval specifying the discount factor. A higher value of delta results in a higher smoothness of the fitted curve. By default, delta is '0.65'.
#' @param prior A list with the prior mean vector \eqn{(m_0)} and covariance matrix \eqn{(C_0)} of \eqn{\theta_0} (state vector at time (age) t = 0). By default mean of zeros and diagonal matrix with a common variance 100 is used. Each element of the list must be named accordingly with the parameter (m0 for mean vector and C0 for covariance matrix).
#' @param prior.sig2 A list with the prior parameters (a, b) of Inverted Gamma distribution for \eqn{\sigma^2}. Each element of the list must be named accordingly with the parameter (a for shape parameter and b for scale parameter).
#' @param M Positive integer indicating the number of iterations of the MCMC run.
#' @param bn Non-negative integer indicating the number of iteration to be discarded as the burn-in period.
#' @param thin Positive integer specifying the period for saving samples.
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
#' random errors with \eqn{v_t \sim N(0, V_t)} and \eqn{w_t \sim N(0, W_t)}. Assumes
#' \eqn{V_t = \sigma^2 V^{'}}, constant for all ages \eqn{t} and \eqn{V^{'}} known. We use
#' the discount factor \eqn{\delta} to specify \eqn{W_t} as \eqn{W_t = C_T(1-\delta)/\delta},
#' where \eqn{C_t} is the conditional covariance matrix of \eqn{\theta_t}. So, if \eqn{\delta = 0}
#' there is no loss information as \eqn{t} increase (completely reducing the smoothness of
#' the fitted curve).
#'
#' To adjust the model adopted a scheme described by (Petris et al, 2009) that uses the FFBS
#' algorithm and a Gibbs sampler step to sample from the parameters of the posterior distribution.
#' For more details, see (Petris et al, 2009).
#'
#' @return A DLM class object.
#' \item{mu}{Posterior samples from \eqn{\mu_t = F_t \theta_t}, for all t.}
#' \item{theta}{Posterior samples from \eqn{\theta_t}, for all t.}
#' \item{sig2}{Posterior samples from \eqn{\sigma^2}.}
#' \item{Wt}{Posterior samples from matrices \eqn{W_t}.}
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
#' qx_t = Dx/Ex
#' qx_t = 1 - exp(-qx_t)
#' y = log(qx_t)
#'
#' ## Fitting DLM
#' fit = dlm(y, M = 100, bn = 20, thin = 1)
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
#' ## credible intervals (See "?qx_ci" for more options):
#' qx_ci(fit)
#'
#' @include ffbs.R
#'
#'@seealso [fitted.DLM()], [predict.DLM()], [plot.DLM()], [print.DLM()] and [summary.DLM()] for `DLM` methods to native R functions [fitted()],
#'[plot()], [print()] and [summary()].
#'
#'[expectancy.DLM()] and [Heatmap.DLM()] for `DLM` methods to compute and visualise the truncated life expectancy
#'via [expectancy()] and [Heatmap()] functions.
#'
#'[dlm_close()] for close methods to expand the life tables.
#'
#'[qx_ci()] and [plot_chain()] to compute credible intervals and visualise the markov chains, respectively.
#'
#' @export
dlm <- function(y, Ft = matrix(c(1,0), nrow = 1), Gt = matrix(c(1,0,1,1), 2), delta = 0.85,
                prior = list(m0 = rep(0, nrow(Gt)), C0 = diag(100, nrow(Gt))),
                prior.sig2 = list(a = 0.01, b = 0.01), M = 5000, bn = 3000, thin = 1,
                ages = 0:(length(y)-1)){

  ## Verificacoes
  if(is.vector(Ft)) {Ft = t(as.matrix(Ft))}
  if(nrow(Ft) != 1) stop("Ft must be a matrix with the following dimensions: 1 row and p columns.")
  if(!(is.matrix(Gt))) {Gt = as.matrix(Gt)}
  if(ncol(Ft) != nrow(Gt)) stop("Matrices Ft and Gt are not well defined.")
  if(ncol(Gt) != nrow(Gt)) stop("Gt must be a square matrix.")
  if(length(prior$m0) != nrow(Gt)) stop("Dimension of prior mean does not match the dimension of matrix Gt.")
  if(nrow(prior$C) != nrow(Gt)) stop("Dimension of prior covariance matrix does not match the dimension of matrix Gt.")
  if(ncol(prior$C) != nrow(Gt)) stop("Dimension of prior covariance matrix does not match the dimension of matrix Gt.")
  if(delta <= 0 || delta >= 1) stop("delta must be in interval (0,1).")

  # Ft0 = matrix(rep(Ft, each = length(y)), ncol = length(Ft))

  ## Initial value for sigma2
  sig2k <- runif(1)

  fit = gibbsSigma2(m0 = prior$m0, C0 = prior$C0, y = as.matrix(y), Ft0 = Ft, Gt = Gt, delta = delta,
                    sig2k = sig2k, alpha = prior.sig2$a, beta = prior.sig2$b, nit = M)

  fit$mu = fit$mu[seq(bn+1, M, by = thin), ]
  fit$theta = fit$theta[seq(bn+1, M, by = thin),, ]
  fit$sig2 = fit$sig2[seq(bn+1, M, by = thin)]
  fit$Wt = fit$Wt[,,,seq(bn+1, M, by = thin)]
  fit$info = list(y = y,
                  ages = ages,
                  Ft = Ft,
                  Gt = Gt,
                  delta = delta,
                  prior = prior,
                  prior.sig2 = prior.sig2)

  return(structure(fit, class = "DLM"))
}
