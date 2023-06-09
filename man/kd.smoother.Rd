% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/kd_smoother.R
\name{kd.smoother}
\alias{kd.smoother}
\title{Kalman Smoother with Drift}
\usage{
kd.smoother(y, f, W, Gt, a2)
}
\arguments{
\item{y}{Observed data as an m-by-N matrix where m is the dimension}

\item{f}{Kalman filtering results as returned by \code{k.filter} on the data}

\item{W}{Covariance p-by-p matrix W referring to system error}

\item{Gt}{Evolution matrix as a p-by-p matrix}

\item{a2}{Drift of the system equation}
}
\value{
A list containing the means \code{s} and the covariances \code{S}
}
\description{
Computes mean and variance of the distribution of the state, conditional
on the covariances of observation and system errors, equation matrices
and all of the observations.
}
\note{
Reference: Petris et al, 2009, p.61
}
