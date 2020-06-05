##' @title Print newline.
##' @param n Integer. Number of lines to print.
##' @return NULL.
##' @author Stephen R. Martin
##' @keywords internal
.newline <- function(n = 1) {
    for(i in 1:n) {
        cat("\n")
    }
}

##' @title Print separator.
##' @param sep String (Default = "-----"). Separator.
##' @param n Integer (Optional). Number of newlines to call after sep.
##' @return NULL.
##' @author Stephen R. Martin
##' @keywords internal
.sep <- function(sep = "\n-----", n = 1) {
    cat(sep)
    if(n > 0) {
        .newline(n)
    }
}

.tab <- function() {
    cat("\t")
}
##' Print method for lmmelsm objects.
##' @title Print method for lmmelsm objects.
##' @param x lmmelsm object.
##' @param ... Not used.
##' @return NULL.
##' @author Stephen R. Martin
##' @export
print.lmmelsm <- function(x, ...) {
    .sep()

    cat("Estimated on: ")
    cat(x$fit@date)
    .newline()
    cat("Elapsed time (seconds): ")
    cat(.get_elapsed_time(x))

    .sep()

    cat("Model specification:")
    .newline(2)

    cat("Measurement model:")
    .newline()
    lapply(x$meta$indicator_spec$mlist, function(f) {
        cat(.print.formula(f))
        .newline()
    })

    if(x$meta$pred_spec$P > 0) {
        .newline()
        cat("Location model:")
        .newline()
        cat(.print.formula(x$meta$pred_spec$plist$location))
        .newline()
    }

    if(x$meta$pred_spec$Q > 0) {
        .newline()
        cat("Scale model")
        .newline()
        cat(.print.formula(x$meta$pred_spec$plist$scale))
        .newline()
    }

    if(x$meta$pred_spec$R > 0) {
        .newline()
        cat("Between-group scale model")
        .newline()
        cat(.print.formula(x$meta$pred_spec$plist$between))
        .newline()
    }
}

summary.lmmelsm <- function(object, prob = .95, ...) {

    dots <- list(...)
    # Define basic structure.
    out <- list(meta = object$meta, summary = list())
    out$meta$digits <- dots$digits %IfNull% 3

    ind_names <- out$meta$indicator_spec$mname
    fnames <- unlist(out$meta$indicator_spec$fname)

    # TODO: Get diagnostics (Rhats, divergences)
    # TODO: Consider: Restructure these as [item/predictor, cols, factor], and name dimensions.
    ## Similar to omegad; makes it easier to print by factor, and to subset.
    ## Or make them 'tidy'

    # Measurement model.
    out$summary$lambda <- .summarize(object, pars = "lambda", prob = prob)
    out$summary$sigma <- .summarize(object, pars = "sigma", prob = prob)
    out$summary$nu <- .summarize(object, pars = "nu", prob = prob)

    ## Restructure
    out$summary$lambda <- .tidy_summary(out$summary$lambda, c("factor", "item"), fnames, ind_names)

    # RE SDs.
    out$summary$mu_logsd_betas_random_sigma <- .summarize(object, pars = "mu_logsd_betas_random_sigma", prob = prob)
    out$summary$Omega_mean_logsd <- .summarize(object, pars = "Omega_mean_logsd", prob = prob)

    # Location predictors
    out$summary$mu_beta <- .summarize(object, pars = "mu_beta", prob = prob)

    # Scale predictors
    out$summary$logsd_beta <- .summarize(object, pars = "logsd_beta", prob = prob)

    # L2 scale predictors
    out$summary$zeta <- .summarize(object, pars = "zeta", prob = prob)

    # Ranefs
    out$summary$mu_random <- .summarize(object, pars = "mu_random", prob = prob)
    out$summary$logsd_random <- .summarize(object, pars = "logsd_random", prob = prob)
    out$summary$mu_beta_random <- .summarize(object, pars = "mu_beta_random", prob = prob)
    out$summary$logsd_beta_random <- .summarize(object, pars = "logsd_beta_random", prob = prob)

    # Eta correlations
    out$summary$Omega_eta <- .summarize(object, pars = "Omega_eta", prob = prob)

    class(out) <- "summary.lmmelsm"
    return(out)
}
##' @title Print method for summary.lmmelsm objects.
##' @param x summary.lmmelsm object.
##' @param ... Not used.
##' @return 
##' @author Stephen Martin
print.summary.lmmelsm <- function(x, ...) {
    dots <- list(...)
    digits <- dots$digits %IfNull% x$meta$digits
}

.print.formula <- function(f) {
    string <- Reduce(paste0, deparse(f))
    return(string)
}
##' @title Gets elapsed time.
##' @param sOut lmmelsm object.
##' @return Numeric.
##' @author Stephen R. Martin
##' @importFrom rstan get_elapsed_time
##' @keywords internal
.get_elapsed_time <- function(sOut) {
    times <- rstan::get_elapsed_time(sOut$fit)
    out <- max(rowSums(times))

    return(out)
}
##' Computes posterior summaries.
##'
##' @title Compute posterior summaries.
##' @param lmmelsm lmmelsm object.
##' @param pars Char vector. Which stan param to summarize.
##' @param prob Numeric (Default: .95; [0 - 1]). The desired mass to contain within the CrI.
##' @return Matrix.
##' @author Stephen R. Martin
##' @keywords internal
.summarize <- function(lmmelsm, pars, prob = .95) {
    samps <- as.matrix(lmmelsm$fit, pars = pars)
    probs <- .prob_to_interval(prob)
    fun <- function(col) {
        m <- mean(col)
        mdn <- quantile(col, probs = .5)
        sd <- sd(col)
        ci <- quantile(col, probs = probs)
        names(ci) <- paste0("Q", probs * 100)
        names(mdn) <- "Median"
        out <- c(Mean = m, mdn, SD = sd, ci)
        return(out)
    }
    samps.sum <- t(apply(samps, 2, fun))
    samps.sum <- as.data.frame(samps.sum)

    return(samps.sum)
}
##' @title Convert stan par-string to numeric columns.
##' @param x String. E.g., "lambda[1,2]"
##' @param labs Character vector (Optional). If supplied, provides the colnames for the matrix.
##' @return Numeric matrix.
##' @author Stephen R. Martin
##' @keywords internal
.pars_to_indices <- function(x, labs = NULL) {
    # inner brackets
    ib <- gsub(".*\\[(.*)\\]", replacement = "\\1", x)
    sep <- strsplit(ib, split = ",")
    num <- lapply(sep, as.numeric)
    inds <- do.call(rbind, num)

    if(!is.null(labs)) { # Apply labels to dimensions.
        colnames(inds) <- labs
    } else { # Label them row, col, arr_1 ... arr_10
        colnames(inds) <- c("row", "col", paste0("arr_", 1:(ncol(inds) - 2)))[1:ncol(inds)]
    }
    
    return(inds)
}

.tidy_summary <- function(x, labs = NULL, ...) {
    dots <- list(...)
    n_relabel <- length(dots)
    
    inds <- .pars_to_indices(rownames(x), labs = labs)
    inds <- as.data.frame(inds)
    for(i in 1:n_relabel) {
        inds[,i] <- dots[[i]][inds[,i]]
    }

    out <- (cbind(inds, x))
    return(out)
}
