% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bmixlm.R
\name{plot.bmixlm}
\alias{plot.bmixlm}
\title{Trace plots for \code{bmixlm} objects}
\usage{
\method{plot}{bmixlm}(x, which = c("probit", "comp1", "comp2", "error"), main = which, ...)
}
\arguments{
\item{x}{An object of class \code{bmixlm}}

\item{which}{The coefficient set to plot}

\item{main}{The main title for the plot}

\item{...}{Additional options to \code{plot.ts}.}
}
\description{
Plot a trace plot of a set of parameters from a \code{bmixlm} object.
}
\details{
The fitted object contains of four sets of parameters: the
coefficients for the two component models, the coefficients for
the probit model, and the error standard deviations from the two
component models. The \code{which} argument determines which of
these parameter sets is plotted:
\itemize{
\item \strong{"comp1"}: coefficients for the first component model
\item \strong{"comp2"}: coefficients for the second component model
\item \strong{"probit"}: coefficients for the probit binary model
\item \strong{"error"}: error standard deviations for the two components
}
}
