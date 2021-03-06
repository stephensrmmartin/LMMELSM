#' Simulated data for fitting the LMMELSM
#'
#' Dataset containing 50 observations of 12 items for 100 persons.
#' The data are generated from an LMMELSM.
#'
#' @format Data frame with 5000 rows and 16 variables.
#' \describe{
#' \item{subject}{The subject ID}
#' \item{baseline}{A subject-level covariate}
#' \item{x1}{A time-varying covariate}
#' \item{x2}{A time-varying covariate}
#' \item{A_1 - A_6}{Six indicators for "Agreeableness"}
#' \item{N_1 - N_6}{Six indicators for "Neuroticism"}
#' }
"sim_data"
