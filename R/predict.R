predict.lmmelsm <- function(object, newdata = NULL, ...) {
    x <- object

    newdata <- newdata %IfNull% x$data

    
}

##' Extract model fitted variates.
##'
##' Extracts model fitted variates. When a latent MMELSM, these are the latent score expectations and log standard deviations. When an observed MMELSM, these are the observed score expectations and log standard deviations.
##' @title Extracted model fitted variates.
##' @param object lmmelsm object.
##' @inheritParams ranef.lmmelsm
##' @param ... Not used.
##' @return List of eta and eta_logsd. If summarize is \code{TRUE}, then the list contains formatted summary tables. If \code{FALSE}, then the list contains MCMC samples for all variables.
##' @author Stephen Martin
fitted.lmmelsm <- function(object, prob = .95, summarize = TRUE, ...) {
    pars <- c("eta", "eta_logsd")
    facVar <- ifelse(object$meta$latent, "factor", "variable")
    fnames <- unlist(object$meta$indicator_spec$fname)

    eta <- as.matrix(object, pars = pars[1])
    eta_logsd <- as.matrix(object, pars = pars[2])

    if(summarize) {
        eta <- .summarize(eta, prob = prob)
        eta_logsd <- .summarize(eta_logsd, prob = prob)

        eta <- .tidy_summary(eta, c("observation", facVar), seq_len(nrow(object$data)), fnames)
        eta_logsd <- .tidy_summary(eta_logsd, c("observation", facVar), seq_len(nrow(object$data)), fnames)

        eta <- .summary_rearrange(eta, c("observation", facVar))
        eta_logsd <- .summary_rearrange(eta_logsd, c("observation", facVar))

        return(nlist(eta, eta_logsd))
    } else {
        return(nlist(eta, eta_logsd))
    }
}
