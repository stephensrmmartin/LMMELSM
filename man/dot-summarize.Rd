% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/print.R
\name{.summarize}
\alias{.summarize}
\title{Compute posterior summaries.}
\usage{
.summarize(lmmelsm, pars, prob = 0.95)
}
\arguments{
\item{lmmelsm}{lmmelsm object.}

\item{pars}{Char vector. Which stan param to summarize.}

\item{prob}{Numeric (Default: .95; \verb{[0 - 1]}). The desired mass to contain within the CrI.}
}
\value{
Matrix.
}
\description{
Computes posterior summaries.
}
\author{
Stephen R. Martin
}
\keyword{internal}
