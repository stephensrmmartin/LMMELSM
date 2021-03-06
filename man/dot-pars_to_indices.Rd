% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/print.R
\name{.pars_to_indices}
\alias{.pars_to_indices}
\title{Convert stan par-string to numeric columns.}
\usage{
.pars_to_indices(x, labs = NULL)
}
\arguments{
\item{x}{String. E.g., \code{"lambda[1,2]"}}

\item{labs}{Character vector (Optional). If supplied, provides the colnames for the matrix.}
}
\value{
Numeric matrix.
}
\description{
Convert stan par-string to numeric columns.
}
\author{
Stephen R. Martin
}
\keyword{internal}
