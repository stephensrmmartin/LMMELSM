% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ranef.R
\name{coef.lmmelsm}
\alias{coef.lmmelsm}
\alias{coef}
\title{Extract group-specific coefficients.}
\usage{
\method{coef}{lmmelsm}(object, prob = 0.95, summarize = TRUE, ...)
}
\arguments{
\item{object}{lmmelsm object.}

\item{prob}{Numeric (Default: .95). Amount of probability mass contained in the credible interval.}

\item{summarize}{Logical (Default: TRUE). Whether to return posterior summaries (TRUE) or MCMC samples (FALSE).}

\item{...}{Not used.}
}
\value{
List of summaries (if \code{summarize} is TRUE), or list of MCMC samples.
}
\description{
Coef method for lmmelsm objects.
}
\details{
Extracts all group-specific coefficients from lmmelsm object.
Note that this is different from \code{\link{ranef}}.
Whereas \code{ranef} extracts the zero-centered random effects, \code{coef} extracts the group-specific effects, defined as the sum of the fixed effect and random effect.
}
\author{
Stephen R Martin
}
