% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_blc_pred.R
\name{plot.PredBLC}
\alias{plot.PredBLC}
\title{BLC: Plot the log-mortality of a prediction}
\usage{
\method{plot}{PredBLC}(x, h = NULL, prob = 0.95, plotIC = TRUE, age = NULL, ...)
}
\arguments{
\item{x}{A \code{PredBLC} object, result to the pred() function call on a \code{BLC} object.}

\item{h}{A numeric vector specifying the year(s) in the prediction horizon to be calculated.}

\item{prob}{A real number that represents the probability of the predict interval.}

\item{plotIC}{Logical. If 'TRUE' (default), shows the predictive intervals.}

\item{age}{A numeric vector indicating the modelled ages. (Optional).}

\item{...}{Other arguments.}
}
\value{
A 'ggplot' object with the predicted mortality rates and their predict intervals.
}
\description{
This functions plot the predicted log-mortality and the predict
intervals of the log-mortality for a specific year in the prediction horizon
}
\examples{
## Importing log-mortality data from Portugal:
data(PT)
Y <- PT

## Fitting the model
fit = blc(Y = Y, M = 100, bn = 20)

#' ## Prediction for 10 years ahead
pred = predict(fit, h = 3)

## Plotting
plot(pred, h = 1)
\donttest{plot(pred, h = 3, prob = 0.9)}

}
\seealso{
\code{\link[=plot.HP]{plot.HP()}}, \code{\link[=plot.DLM]{plot.DLM()}} and \link{plot.BLC} for \code{HP}, \code{DLM} or \code{BLC} methods.
}
