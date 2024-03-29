% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hp_close.R
\name{hp_close}
\alias{hp_close}
\title{HP: Fitting the advanced ages of the life tables}
\usage{
hp_close (fit, method = c("hp", "plateau", "linear", "gompertz"),
 x0 = max(fit$data$x), max_age = 120, k = 7,
 weights = seq(from = 0, to = 1, length.out = 2*k+1),
 new_Ex = NULL, new_Dx = NULL)
}
\arguments{
\item{fit}{Object of the class \code{HP} fitted by the hp() function}

\item{method}{Character string specifying the closing method to be fitted, with them being: 'hp', 'plateau', 'linear' or 'gompertz'.}

\item{x0}{Integer with the starting age the closing method will be fitted from. Default is the last age fitted by the 'HP' object.}

\item{max_age}{Integer with the maximum age the closing method will be fitted. Default age is '120'.}

\item{k}{Integer representing the size of the age interval to be mixed with the 'linear' or 'gompertz' closing methods for smooth graduation. If k = 0, no mixing will be applied.}

\item{weights}{Vector of weights of the closing method used in the mixture of the closing method and the fitted model made in the mixing age group. The vector's size should be equal to 2k+1. For a better understanding of this parameter and the mixture applied in this function, see Details.}

\item{new_Ex}{Vector with exposure of ages after the x0 input. This is an optional argument used in the 'linear' and 'gompertz' closing methods. If this argument is specified, then new_Dx also needs to be.}

\item{new_Dx}{Vector containing the death counts of the ages after the x0 input. This is also an optional argument used in the 'linear' and 'gompertz' closing methods. The length must be the same as new_Ex.}
}
\value{
Returns a \code{ClosedHP} class object with the predictive chains of the death probability
(qx) from first fitted age to max_age argument, the data utilized by the function and the
closing method chosen.
}
\description{
This function receives an object of the class \code{HP} fitted by the hp() function
and fits a closing method to expand the life tables dataset to a maximum age argument inputted
by the user.
There are four closing methods available: 'hp', 'plateau', 'linear', and 'gompertz'.
The 'linear' method can only be applied with HP objects following the lognormal variant of
the HP mortality law.
}
\details{
There are three types of age groups when the closing method is applied: a group
where only the HP-fitted model computes the death probabilities, followed by a
group in which the death probabilities are a mix (or more precise a weighted mean)
from the HP model and the closing method and followed by a group in which the
death probabilities are computed just by the closing method. The mix is applied
so the transition of the death probabilities of the ages between the fitted model
and the closing method occurs smoothly.

The parameters 'x0' and 'k' define the mixing group age. The parameter 'x0'
indicates the center age of the group. The parameter 'k' is the range of ages
before 'x0' and after 'x0', so this group has a total of \eqn{2k + 1} age. Therefore,
the parameter 'weights' must have a length size equal to \eqn{2k + 1}. In this case,
the death probability is calculated as follows. Consider \eqn{model_x} and \eqn{close_x}
as the death probability of the fitted model and closing method in the age \eqn{x},
respectively. Then, the resulting death probability of the mix is calculated as:

\eqn{q_x = w_x model_x + (1-w_x)close_x},

where \eqn{w_x} represents the weight of the closing method in the age \eqn{x}.
This computation is applied to all elements in the MCMC chain of the fitted model,
resulting in a new chain of death probabilities. This procedure is applied only in
the linear and Gompertz methods.

The four closing methods for life tables are:

1.'hp' method: Expands the previously adjusted HP model until the max_age argument.

2.'plateau' method: Keeps the death probability (qx) constant after the x0 argument.

3.'linear' method: Fits a linear regression starting at age x0 - k until the last age with data available (lognormal only).

4.'gompertz' method: Adopted as the closing method of the 2010-2012 English Life Table No. 17, fits the Gompertz mortality law via SIR using the same available data as the 'linear' method.
}
\examples{
## Importing mortality data from the USA available on the Human Mortality Database (HMD):
data(USA)

## Selecting the exposure and the death count of the year 2010, ranging from 0 to 90 years old:
USA2010 = USA[USA$Year == 2010,]
x = 0:90
Ex = USA2010$Ex.Male[x+1]
Dx = USA2010$Dx.Male[x+1]

## Fitting a lognormal HP model:
fit = hp(x = x, Ex = Ex, Dx = Dx, model = "lognormal",
         M = 1000, bn = 0, thin = 10)

## Applying the closing function with different methods:
close1 = hp_close(fit, method = "hp", x0 = 90)
\donttest{close2 = hp_close(fit, method = "plateau", x0 = 90)
close3 = hp_close(fit, method = "linear", x0 = 80,
                  new_Ex = USA2010$Ex.Male[82:101],
                  new_Dx = USA2010$Dx.Male[82:101])
close4 = hp_close(fit, method = "gompertz", x0 = 70,
                  new_Ex = USA2010$Ex.Male[72:101],
                  new_Dx = USA2010$Dx.Male[72:101],
                  k = 5, max_age = 120)

#### Using the other functions available in the package with the 'ClosedHP' object:

## qx estimation (See "?fitted.HP" in the BayesMortalityPlus package for more options):
fitted(close2)

## life expectancy (See "?expectancy.HP" for more options)
expectancy(close3, age = 0:110)

## plotting (See "?plot.HP" in the BayesMortalityPlus package for more options):
plot(close4)
g <- plot(list(close4, fit),
          colors = c("seagreen", "blue"),
          labels = c("Closed", "Model"))
# plotly::ggplotly(g)
}

}
\references{
Dodd, Erengul, Forster, Jonathan, Bijak, Jakub, & Smith, Peter 2018. “Smoothing mortality data: the English life table, 2010-12.” \emph{Journal of the Royal Statistical Society: Series A (Statistics in Society)}, 181(3), 717-735.
}
\seealso{
\code{\link[=fitted.HP]{fitted.HP()}}, \code{\link[=plot.HP]{plot.HP()}}, \code{\link[=print.HP]{print.HP()}} and \code{\link[=summary.HP]{summary.HP()}} for \code{ClosedHP} methods to native R functions \code{\link[=fitted]{fitted()}},
\code{\link[=plot]{plot()}}, \code{\link[=print]{print()}} and \code{\link[=summary]{summary()}}.

\code{\link[=expectancy.HP]{expectancy.HP()}} and \code{\link[=Heatmap.HP]{Heatmap.HP()}} for \code{ClosedHP} methods to compute and visualise the truncated life expectancy
via \code{\link[=expectancy]{expectancy()}} and \code{\link[=Heatmap]{Heatmap()}} functions.
}
