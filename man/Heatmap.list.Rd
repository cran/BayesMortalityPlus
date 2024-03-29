% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/heatmap_list.R
\name{Heatmap.list}
\alias{Heatmap.list}
\title{Heatmap for a set of life tables}
\usage{
\method{Heatmap}{list}(
  x,
  x_lab = NULL,
  age = NULL,
  max_age = NULL,
  color = c("red", "white", "blue"),
  ...
)
}
\arguments{
\item{x}{List of objects of classes: \code{HP}, \code{DLM}, \code{ClosedHP}, or \code{ClosedDLM}.}

\item{x_lab}{Description of the object 'fit'.}

\item{age}{Vector with the ages to plot the heatmap.}

\item{max_age}{Positive number indicating the last age to be considered to compute the life expectancy (methods for matching the age interval will be considered if needed). This argument is only necessary with objects of the class \code{HP} or \code{DLM}.}

\item{color}{Vector with colours used in the heatmap.}

\item{...}{Further arguments passed to or from other methods.}
}
\value{
A ggplot2 heatmap of the life expectancy.
}
\description{
This function plots a heatmap for the life expectancy of the mortality graduations
returned by hp(), dlm(), hp_close() or dlm_close() functions.
}
\examples{
## Importing mortality data from the USA available on the Human Mortality Database (HMD):
data(USA)

# Example (HP): -------------------------------

## Selecting the data from 2010
USA2010 = USA[USA$Year == 2010,]

ExF = USA2010$Ex.Female[1:91]
DxF = USA2010$Dx.Female[1:91]
x <- 0:90

fitF <- hp(x, ExF, DxF, model = "lognormal", M = 1000, bn = 0, thin = 10)

ExM = USA2010$Ex.Male[1:91]
DxM = USA2010$Dx.Male[1:91]

fitM <- hp(x, ExM, DxM, model = "lognormal", M = 1000, bn = 0, thin = 10)

fits <- list(fitF = fitF, fitM = fitM)

Heatmap(fits, x_lab = c("Female 2010 USA","Male 2010 USA"),
        age = 15:85)


}
\seealso{
\code{\link[=Heatmap.HP]{Heatmap.HP()}}, \code{\link[=Heatmap.DLM]{Heatmap.DLM()}} and \code{\link[=Heatmap.BLC]{Heatmap.BLC()}} for drawing single Heatmaps.
}
