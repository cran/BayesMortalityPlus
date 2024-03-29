% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/print_dlm.R
\name{print.DLM}
\alias{print.DLM}
\title{DLM: Print}
\usage{
\method{print}{DLM}(x, ...)
}
\arguments{
\item{x}{A \code{DLM} or \code{ClosedDLM} object, result of a call to dlm() or dlm_close() function.}

\item{...}{Further arguments passed to or from other methods.}
}
\value{
A character vector with the details of a fitted \code{DLM} or \code{ClosedDLM} model.
}
\description{
Print details from a fitted \code{DLM} or \code{ClosedDLM} models and returns it invisibly.
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
fit = dlm(y, M = 100)
print(fit)

}
\seealso{
\code{\link[=print.HP]{print.HP()}} and \code{\link[=print.BLC]{print.BLC()}} for \code{HP} or \code{BLC} methods.
}
