% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mean_blc_pred.R
\name{mean.PredBLC}
\alias{mean.PredBLC}
\title{BLC: Arithmetic mean for predictions}
\usage{
\method{mean}{PredBLC}(x, h, ...)
}
\arguments{
\item{x}{A \code{PredBLC} object, result to the pred() function call on a \code{BLC} object.}

\item{h}{A positive integer specifying the year in the prediction horizon to be calculated.}

\item{...}{Further arguments passed to or from other methods.}
}
\value{
A vector with the mean values of the log-mortality chains.
}
\description{
Calculates the means based on the resulting chains from a predicted year.
}
\examples{
data(PT)
Y <- PT

## Fitting the model
fit = blc(Y = Y, M = 100, bn = 20)

## Prediction for 2 years ahead
pred = predict(fit, h = 2)

mean(pred, 1)
mean(pred, 2)

}
\seealso{
\code{\link[=mean.BLC]{mean.BLC()}} for \code{BLC} object method.
}
