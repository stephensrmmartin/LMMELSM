
##' @title loo method for LMMELSM objects.
##' @param x lmmelsm object.
##' @param ... Not used.
##' @return loo object.
##' @author Stephen R. Martin
##' @importFrom loo loo
##' @export loo
##' @export
loo.lmmelsm <- function(x, ...) {
    
}

.log_liks <- function(x) {
    # Extract and transform param samples.
    etas <- .extract_transform(x$fit, par = "eta")
    lambdas <- .extract_transform(x$fit, par = "lambda")
    nus <- .extract_transform(x$fit, par = "nu")
    sigmas <- .extract_transform(x$fit, par = "sigma")

    
}
