
##' @title loo method for LMMELSM objects.
##' @param x lmmelsm object.
##' @param type String (Default: "observation"). If "observation", then loo is leave-row-out. If "group", then loo is leave-group-out.
##' @param ... Not used.
##' @return loo object.
##' @author Stephen R. Martin
##' @importFrom loo loo
##' @export loo
##' @aliases loo
##' @export
loo.lmmelsm <- function(x, type = c("observation", "group"), ...) {
    type <- match.arg(type)

    ll <- .log_liks(x)

    ll_fun <- switch(type,
                     observation = .ll_obs,
                     group = .ll_group)

    ll <- ll_fun(ll, x)

    chain_id <- rep(1:x$stan_args$chains, each = x$stan_args$iter - x$fit@stan_args[[1]]$warmup)
    r_eff <- loo::relative_eff(exp(ll), chain_id = chain_id)
    looOut <- loo::loo(ll, r_eff = r_eff)

    return(looOut)
}

.ll_obs <- function(ll, ...) {

    out <- apply(ll, 3, function(x) { # Each sample
        rowSums(x) # Observation-LL; 
    }) # N x S
    out <- t(out) # Make S x N

    return(out)
    
}

.ll_group <- function(ll, x) {
    ll_obs <- .ll_obs(ll)
    ll_obs <- t(ll_obs) # Make N x S

    gn <- x$meta$group_spec$numeric
    N <- nrow(ll_obs)
    N_seq <- seq_len(N)
    ## out <- tapply(ll_obs, gn, colSums)
    out <- tapply(N_seq, gn, function(x) {
        colSums(ll_obs[x,])
    }, simplify = FALSE) # K-length list; ordered by 1:K
    out <- do.call(rbind, out) # K x S
    out <- t(out)  # S x K

    return(out)
}

.log_liks <- function(x) {
    J <- get_J(x)
    N <- get_N(x)
    y <- x$stan_data$y
    # Extract and transform param samples.
    etas <- .extract_transform(x$fit, par = "eta")
    lambdas <- .extract_transform(x$fit, par = "lambda")
    nus <- .extract_transform(x$fit, par = "nu")
    sigmas <- .extract_transform(x$fit, par = "sigma")

    S <- dim(etas)[3]

    # Predictions
    preds <- array(0, dim = c(N, J, S))
    for(s in seq_len(S)) {
        ## preds[,, s] <- .ones(N) %*% t(nus[, s]) + etas[,, s] %*% lambdas[,, s]
        preds[,, s] <- .ones(N) %*% t(nus[, s]) + .array_extract(etas, s) %*% .array_extract(lambdas, s)
    }

    # Likelihoods
    log_lik <- array(0, dim = c(N, J, S))
    for(n in seq_len(N)) {
        log_lik[n, ,] <- dnorm(y[n, ], preds[n, , ], sigmas, log = TRUE)
    }

    return(log_lik)
}
