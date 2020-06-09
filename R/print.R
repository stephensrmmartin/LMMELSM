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
    out$summary[c("lambda", "sigma", "nu")] <- .summary_measurement(object, prob)

    # Random effects
    out$summary[c("mu_logsd_betas_random_sigma",
                  "Omega_mean_logsd",
                  "mu_random",
                  "logsd_random",
                  "mu_beta_random",
                  "logsd_beta_random")] <- .summary_ranef(object, prob)

    # Fixed effects
    out$summary[c("mu_beta", "logsd_beta")] <- .summary_fixef(object, prob)

    # L2 scale predictors
    out$summary$zeta <- .summary_between(object, prob)

    # Eta correlations
    out$summary$Omega_eta <- .summary_epsilon(object, prob)

    class(out) <- "summary.lmmelsm"
    return(out)
}

##' @title Rearrange summary output.
##' @param x Summary table.
##' @param cols Columns, in order, to place in front.
##' @param arrange Logical (Default: TRUE). Whether to sort rows.
##' @return Data.frame.
##' @author Stephen R. Martin
##' @keywords internal
.summary_rearrange <- function(x, cols, arrange = TRUE) {
    cns <- colnames(x)
    col_ind <- seq_len(ncol(x))
    where <- match(cols, cns)
    newcols <- c(where, col_ind[-where])
    x <- x[, newcols]

    if(arrange) {
        args <- x[, cols]
        ord <- do.call(order, args)
        x <- x[ord, ]
    }

    return(x)
}
# TODO: Fix the .summary_* fns using .summary_rearrange.

.summary_measurement <- function(x, prob) {
    fnames <- unlist(x$meta$indicator_spec$fname)
    ind_names <- x$meta$indicator_spec$mname

    lambda <- .summarize(x, pars = "lambda", prob = prob)
    sigma <- .summarize(x, pars = "sigma", prob = prob)
    nu <- .summarize(x, pars = "nu", prob = prob)

    lambda <- .tidy_summary(lambda, c("factor", "item"), fnames, ind_names)
    sigma <- .tidy_summary(sigma, "item", ind_names)
    nu <- .tidy_summary(nu, "item", ind_names)

    out <- nlist(lambda, sigma, nu)
    return(out)
}

.summary_ranef <- function(x, prob) {
    # Meta-data
    pnames <- x$meta$pred_spec$pname
    fnames <- unlist(x$meta$indicator_spec$fname)
    F <- x$meta$indicator_spec$F
    P_random <- x$meta$pred_spec$P_random
    Q_random <- x$meta$pred_spec$Q_random
    P_random_ind <- x$meta$pred_spec$P_random_ind
    Q_random_ind <- x$meta$pred_spec$Q_random_ind

    # Build RE names
    re_names <- paste0(fnames, "MAGICSEP", rep(c("mu", "logsd"), each = F))
    re_total <- 2 * F + 2 * P_random + 2 * Q_random

    if(P_random > 0) {
        re_names <- c(re_names, paste0(rep(fnames, each = P_random), "MAGICSEP", pnames$location[P_random_ind]))
    }
    if(Q_random > 0) {
        re_names <- c(re_names, paste0(rep(fnames, each = Q_random), "MAGICSEP", pnames$scale[Q_random_ind]))
    }
    re_int_names <- re_names[1:(2 * F)]
    re_slope_names <- re_names[(2 * F + 1):re_total]

    ###############
    # Hyperpriors #
    ###############
    # Summarize
    mu_logsd_betas_random_sigma <- .summarize(x, pars = "mu_logsd_betas_random_sigma", prob = prob)
    Omega_mean_logsd <- .summarize(x, pars = "Omega_mean_logsd", prob = prob)

    # Tidy
    mu_logsd_betas_random_sigma <- .tidy_summary(mu_logsd_betas_random_sigma, "param", re_names)
    Omega_mean_logsd <- .tidy_summary(Omega_mean_logsd, c("row", "col"), re_names, re_names)

    # Separate strings
    mu_logsd_betas_random_sigma[, c("factor", "param")] <- .magicsep(mu_logsd_betas_random_sigma[, "param"], c("factor", "param"))

    Omega_mean_logsd[, c("row_factor", "row_param")] <- .magicsep(Omega_mean_logsd[, "row"], c("row_factor", "row_param"))
    Omega_mean_logsd[, c("col_factor", "col_param")] <- .magicsep(Omega_mean_logsd[, "col"], c("col_factor", "col_param"))
    Omega_mean_logsd[, c("row", "col")] <- NULL
    # Only return lower.tri
    Omega_mean_logsd <- Omega_mean_logsd[.full_to_lower_tri(re_total), ]

    ##################
    # Random Effects #
    ##################
    gs <- x$meta$group_spec
    # Summarize
    mu_random <- .summarize(x, pars = "mu_random", prob = prob)
    logsd_random <- .summarize(x, pars = "logsd_random", prob = prob)
    mu_beta_random <- .summarize(x, pars = "mu_beta_random", prob = prob)
    logsd_beta_random <- .summarize(x, pars = "logsd_beta_random", prob = prob)

    # Tidy
    mu_random <- .tidy_summary(mu_random, c(gs$name, "factor"), gs$map$label, fnames)
    logsd_random <- .tidy_summary(logsd_random, c(gs$name, "factor"), gs$map$label, fnames)
    mu_beta_random <- .tidy_summary(mu_beta_random, c(gs$name, "predictor", "factor"), gs$map$label, pnames$location[P_random_ind], fnames)
    logsd_beta_random <- .tidy_summary(logsd_beta_random, c(gs$name, "predictor", "factor"), gs$map$label, pnames$scale[Q_random_ind], fnames)

    out <- nlist(mu_logsd_betas_random_sigma,
                 Omega_mean_logsd,
                 mu_random,
                 logsd_random,
                 mu_beta_random,
                 logsd_beta_random
                 )

    return(out)
}

.summary_fixef <- function(x, prob) {
    # Meta-data
    pnames <- x$meta$pred_spec$pname
    fnames <- unlist(x$meta$indicator_spec$fname)

    # Summarize
    mu_beta <- .summarize(x, pars = "mu_beta", prob = prob)
    logsd_beta <- .summarize(x, pars = "logsd_beta", prob = prob)

    # Tidy
    mu_beta <- .tidy_summary(mu_beta, c("predictor", "factor"), pnames$location, fnames)
    logsd_beta <- .tidy_summary(logsd_beta, c("predictor", "factor"), pnames$scale, fnames)

    out <- nlist(mu_beta, logsd_beta)
    return(out)
}

.summary_between <- function(x, prob) {
    # Meta-data
    pnames <- x$meta$pred_spec$pname$between
    fnames <- unlist(x$meta$indicator_spec$fname)
    F <- x$meta$indicator_spec$F
    re_int_names <- paste0(fnames, "MAGICSEP", rep(c("mu", "logsd"), each = F))

    # Summarize
    zeta <- .summarize(x, pars = "zeta", prob = prob)

    # Tidy
    zeta <- .tidy_summary(zeta, c("predictor", "param"), pnames, re_int_names)

    # Separate
    zeta[, c("factor", "param")] <- .magicsep(zeta[, "param"], c("factor", "param"))

    out <- nlist(zeta)
    return(out)
}

.summary_epsilon <- function(x, prob) {
    # Meta-data
    fnames <- unlist(x$meta$indicator_spec$fname)
    F <- x$meta$indicator_spec$F

    # Summarize
    o <- .summarize(x, pars = "Omega_eta", prob = prob)

    # Tidy
    o <- .tidy_summary(o, c("row", "col"), fnames, fnames)

    # Lower-tri only
    o <- o[.full_to_lower_tri(F), ]

    out <- list(Omega_eta = o)
    return(out)
}
##' @title Print method for summary.lmmelsm objects.
##' @param x summary.lmmelsm object.
##' @param ... Not used.
##' @return 
##' @author Stephen R. Martin
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

##' Creates "tidy" summaries in lieu of the stan rownames.
##'
##' .summarize creates an rstan-like summary with rownames, mat[1:R, 1:C].
##' \code{.tidy_summary(mat, c("rows", "cols"))} would then create two new columns, "rows" and "cols" with the indices in them.
##' If arguments are provided in \code{...}, then these indicate the mappings between the indices and labeled values.
##' E.g., \code{.tidy_summary(mat, c("rows", "cols"), c("A", "B"), c("C", "D"))} would create two new columns, "rows" and "cols", and replace rows = 1 with rows = A; cols=2 with cols = D, and so on.
##' Useful for going from stan rownames, to labeled columns.
##' @title Takes stan summary, returns summary with indices-as-columns.
##' @param x Output of .summarize
##' @param labs The labels for each parameter index. E.g., "predictor", "factor"
##' @param ... Optional (but recommended). Mappings for indices. E.g., Index column 1 is replaced by ...[[1]][col1Indices].
##' @return Data frame.
##' @author Stephen R. Martin
##' @keywords internal
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
##' Converts a character vector int columns.
##'
##' A wrapper around strcapture.
##' Given "a_MAGICSEP_b", returns a vector of "a_" and "_b".
##' Useful for converting a multiple-parameter string into columns of parameters.
##' User can give column names in \code{labs}, and types of each extracted component as a list in \code{types}.
##' User can customize \code{sep}, but by default assumes MAGICSEP.
##' @title Convert char vector to columns.
##' @param charvec Character vector. Characters to separate into columns.
##' @param labs Character vector. Labels for the columns.
##' @param types List (Default: All characters). List (in order) of extracted types.
##' @param sep String (Default: "MAGICSEP").
##' @return data.frame
##' @author Stephen R. Martin
##' @keywords internal
.magicsep <- function(charvec, labs, types = NULL, sep = "MAGICSEP") {
    args <- list()
    if(is.null(types)) {
        for(l in labs) {
            args[[l]] <- character()
        }
    } else {
        for(i in seq_len(length(labs))) {
            args[[labs[i]]] <- types[[i]]
        }
    }
    proto <- do.call(data.frame, args)
    out <- strcapture(paste0("(.*)",sep,"(.*)"), charvec, proto)
    return(out)
}
##' Helper for correlation-matrix summarize output.
##'
##' The .summarize function returns every redundant and constant element from a correlation matrix.
##' This function returns the stan-strings (when \code{string = TRUE}, e.g., "[2,1]", "[3,1]"), or the row-index assuming column-major order.
##' @title Get indices for subsetting lower-tri summaries of square matrices. 
##' @param x Integer. Dimension of matrix.
##' @param string Logical (Default: FALSE). Whether to return strings (e.g., "[2,1]", or row indices, assuming column-major ordering.)
##' @return 
##' @author Stephen R. Martin
##' @keywords internal
.full_to_lower_tri <- function(x, string = FALSE) {
    if(string) {
        inds <- which(lower.tri(matrix(0, x, x)), arr.ind = TRUE)
        out <- paste0("[", inds[, 1], ",", inds[,2], "]")
        return(out)
    } else {
        inds <- which(lower.tri(matrix(0, x, x)))
        out <- inds
        return(out)
    }
}
