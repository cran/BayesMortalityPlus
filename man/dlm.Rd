% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dlm.R
\name{dlm}
\alias{dlm}
\title{Dynamic Linear Model for mortality table graduation}
\usage{
dlm(y, Ft = matrix(c(1,0), nrow = 1), Gt = matrix(c(1,0,1,1), 2), delta = 0.85,
    prior = list(m0 = rep(0, nrow(Gt)), C0 = diag(100, nrow(Gt))),
    prior.sig2 = list(a = 0, b = 0), M = 2000, ages = 0:(length(y)-1))
}
\arguments{
\item{y}{Numeric vector of log mortality rates.}

\item{Ft}{1xp Matrix that specifies the observation equation, where p is the number of parameters. By default, 'Ft = matrix(c(1,0), nrow = 1)'.}

\item{Gt}{pxp Matrix that specifies the system equations. By default, Gt = matrix(c(1,0,1,1), 2).}

\item{delta}{Positive real value or real vector of the same length as y with values in the '(0, 1)' interval specifying the discount factor for each age. A higher value of delta results in a higher smoothness of the fitted curve. If a single value is defined, this same value is used for all ages. By default, delta is '0.85'.}

\item{prior}{A list with the prior mean vector \eqn{(m_0)} and covariance matrix \eqn{(C_0)} of \eqn{\theta_0} (state vector at time (age) t = 0). By default mean of zeros and diagonal matrix with a common variance 100 is used. Each element of the list must be named accordingly with the parameter (m0 for mean vector and C0 for covariance matrix).}

\item{prior.sig2}{A list with the prior parameters (a, b) of Inverted Gamma distribution for \eqn{\sigma^2}. Each element of the list must be named accordingly with the parameter (a for shape parameter and b for scale parameter).}

\item{M}{Positive integer that indicates the sampling size from the posterior distributions. The default value is 2000.}

\item{ages}{Numeric vector of the ages fitted. Default is '0:(length(y)-1)'.}
}
\value{
A DLM class object.
\item{mu}{Posterior samples from \eqn{\mu_t = F_t \theta_t}, for all t.}
\item{theta}{Posterior samples from \eqn{\theta_t}, for all t.}
\item{sig2}{Posterior samples from \eqn{\sigma^2}.}
\item{param}{A list with the states parameters for filtering distribution (mt, Ct), predictive distribution (ft, Qt), smoothing distribution (as, Rs), and parameters of the posterior distribution for variance (alpha, beta).}
\item{info}{A list with some informations of the fitted model: the specification of \eqn{F_t} and \eqn{G_t} matrices, the data y and the ages, the discount factor \eqn{delta} value specified and priors informations.}
}
\description{
This function fits a Dynamic Linear Model (DLM) for mortality data following
a Bayesian framework using Forward Filtering Backward Sampling algorithm to compute the posterior distribution.
The response variable is the log of the mortality rate, and it is modeled specifying the matrices Ft and Gt from the DLM equations.
Furthermore, the discount factor is used to control the smoothness of the fitted model. By default, a
linear growth model is specified.
}
\details{
Let \eqn{Y_t} be the log mortality rate at age \eqn{t}. A DLM is specified as follows:

For \eqn{t = 0}:

\eqn{\theta_0 \sim N_p (m_0, C_0)}

Now, for \eqn{t \geq 1}:

The observation equation:

\eqn{Y_t = F_t \theta_t + v_t}

The system equation:

\eqn{\theta_t = G_t \theta_{t-1} + w_t}

Where \eqn{F_t} and \eqn{G_t} are known matrices. \eqn{v_t} and \eqn{w_t} are independent
random errors with \eqn{v_t \sim N(0, \sigma^2)} and \eqn{w_t \sim N(0, \sigma^2 W_t)}. We
use the discount factors \eqn{\delta} to specify \eqn{W_t} as \eqn{W_t = C_t(1-\delta)/\delta},
where \eqn{C_t} is the conditional covariance matrix of \eqn{\theta_t}. So, if
\eqn{\delta = 0} there is no loss information as \eqn{t} increase (completely reducing the
smoothness of the fitted curve). \eqn{\delta} can be specified as a single value for all ages
or as a vector in which each element is associated with an age.

A scheme described by (Petris et al, 2009) for conjugated inference is used.
For more details, see (Petris et al, 2009).
}
\examples{
## Importing mortality data from the USA available on the Human Mortality Database (HMD):
data(USA)

## Selecting the log mortality rate of the 2010 male population ranging from 0 to 100 years old
USA2010 = USA[USA$Year == 2010,]
x = 0:100
Ex = USA2010$Ex.Male[x+1]
Dx = USA2010$Dx.Male[x+1]
y = log(Dx/Ex)

## Fitting DLM
fit = dlm(y)
print(fit)
summary(fit)

## Using other functions available in the package:
## plotting (See "?plot.DLM" in the BayesMortality package for more options):
plot(fit)

## qx estimation (See "?fitted.DLM" in the BayesMortality package for more options):
fitted(fit)

## chain's plot (See "?plot_chain" for more options):
plot_chain(fit, param = c("mu[0]", "mu[100]"))

## Varying discount factor
fit2 = dlm(y, delta = c(rep(0.8, 36), rep(0.9, 65)))
plot(fit2)


}
\references{
Campagnoli, P., Petris, G., and Petrone, S. (2009). \emph{Dynamic linear models with R}. Springer-Verlag New York.
}
\seealso{
\code{\link[=fitted.DLM]{fitted.DLM()}}, \code{\link[=predict.DLM]{predict.DLM()}}, \code{\link[=plot.DLM]{plot.DLM()}}, \code{\link[=print.DLM]{print.DLM()}} and \code{\link[=summary.DLM]{summary.DLM()}} for \code{DLM} methods to native R functions \code{\link[=fitted]{fitted()}},
\code{\link[=plot]{plot()}}, \code{\link[=print]{print()}} and \code{\link[=summary]{summary()}}.

\code{\link[=expectancy.DLM]{expectancy.DLM()}} and \code{\link[=Heatmap.DLM]{Heatmap.DLM()}} for \code{DLM} methods to compute and visualise the truncated life expectancy
via \code{\link[=expectancy]{expectancy()}} and \code{\link[=Heatmap]{Heatmap()}} functions.

\code{\link[=dlm_close]{dlm_close()}} for close methods to expand the life tables.

\code{\link[=plot_chain]{plot_chain()}} to visualise the markov chains, respectively.
}
