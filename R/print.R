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
##' @return x (Invisibly).
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
    invisible(x)
}
##' @title Summary method for lmmelsm objects.
##' @param object lmmelsm object.
##' @param prob Numeric (Default: .95). Amount of probability mass contained in the credible interval.
##' @param ... Not used.
##' @return summary.lmmelsm object. A list containing \code{meta} (metadata) and \code{summary} (summary tables).
##' @author Stephen R. Martin
##' @method summary lmmelsm
##' @export
summary.lmmelsm <- function(object, prob = .95, ...) {

    dots <- list(...)
    # Define basic structure.
    out <- list(meta = object$meta, summary = list())
    out$meta$digits <- dots$digits %IfNull% 3

    # More meta-data
    out$meta$stan <- list()
    out$meta$stan$date <- object$fit@date
    out$meta$stan$elapsed <- rstan::get_elapsed_time(object$fit)
    out$meta$stan$diag <- .get_diagnostics(object)

    latent <- object$meta$latent
    if(latent) {
        # Measurement model.
        out$summary[c("lambda", "sigma", "nu")] <- .summary_measurement(object, prob)
    } else {
        out$summary[c("sigma", "nu")] <- .summary_observed(object, prob)
    }

    # Random effects
    out$summary[c("random_sigma",
                  "random_correlation",
                  "random_mu_intercept",
                  "random_logsd_intercept",
                  "random_mu_coef",
                  "random_logsd_coef")] <- .summary_ranef(object, prob)

    # Fixed effects
    out$summary[c("mu_coef", "logsd_coef")] <- .summary_fixef(object, prob)

    # L2 scale predictors
    out$summary$zeta <- .summary_between(object, prob)$zeta

    # Eta correlations
    out$summary$factor_correlation <- .summary_epsilon(object, prob)$Omega_eta

    class(out) <- "summary.lmmelsm"
    return(out)
}

##' @title Rearrange summary output.
##' @param x Summary table.
##' @param cols Columns, in order, to place in front.
##' @param arrange Logical (Default: FALSE). Whether to sort rows.
##' @return Data.frame.
##' @author Stephen R. Martin
##' @keywords internal
.summary_rearrange <- function(x, cols, arrange = FALSE) {
    cns <- colnames(x)
    col_ind <- seq_len(ncol(x))
    where <- match(cols, cns)
    newcols <- c(where, col_ind[-where])
    x <- x[, newcols]

    if(arrange) {
        args <- x[, cols, drop = FALSE]
        ord <- do.call(order, args)
        x <- x[ord, ]
    }

    return(x)
}

.summary_measurement <- function(x, prob) {
    # Meta-data
    fnames <- unlist(x$meta$indicator_spec$fname)
    ind_names <- x$meta$indicator_spec$mname

    # Summarize
    lambda <- .summarize(x, pars = "lambda", prob = prob)
    sigma <- .summarize(x, pars = "sigma", prob = prob)
    nu <- .summarize(x, pars = "nu", prob = prob)

    # Tidy
    lambda <- .tidy_summary(lambda, c("factor", "item"), fnames, ind_names)
    sigma <- .tidy_summary(sigma, "item", ind_names)
    nu <- .tidy_summary(nu, "item", ind_names)

    # Rearrange
    lambda <- .summary_rearrange(lambda, c("factor", "item"))
    sigma <- .summary_rearrange(sigma, "item")
    nu <- .summary_rearrange(nu, "item")

    out <- nlist(lambda, sigma, nu)
    return(out)
}

.summary_observed <- function(x, prob) {
    # Meta-data
    ind_names <- x$meta$indicator_spec$mname

    # Summarize
    sigma <- .summarize(x, pars = "sigma", prob = prob)
    nu <- .summarize(x, pars = "nu", prob = prob)

    # Tidy
    sigma <- .tidy_summary(sigma, "item", ind_names)
    nu <- .tidy_summary(nu, "item", ind_names)

    # Rearrange
    sigma <- .summary_rearrange(sigma, "item")
    nu <- .summary_rearrange(nu, "item")

    out <- nlist(sigma, nu)
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
    re_total <- 2 * F + F*P_random + F*Q_random
    re_names <- .build_re_names(x$meta)
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

    # Rearrange
    mu_logsd_betas_random_sigma <- .summary_rearrange(mu_logsd_betas_random_sigma, c("factor", "param"), arrange=FALSE)
    Omega_mean_logsd <- .summary_rearrange(Omega_mean_logsd, c("row_factor", "row_param", "col_factor", "col_param"), arrange=FALSE)

    ##################
    # Random Effects #
    ##################
    gs <- x$meta$group_spec

    mu_random <- .summarize(x, pars = "mu_random", prob = prob)
    mu_random <- .tidy_summary(mu_random, c(gs$name, "factor"), gs$map$label, fnames)
    mu_random <- .summary_rearrange(mu_random, c(gs$name, "factor"))

    logsd_random <- .summarize(x, pars = "logsd_random", prob = prob)
    logsd_random <- .tidy_summary(logsd_random, c(gs$name, "factor"), gs$map$label, fnames)
    logsd_random <- .summary_rearrange(logsd_random, c(gs$name, "factor"))

    mu_beta_random <- logsd_beta_random <- NA

    if(P_random > 0) {
        mu_beta_random <- .summarize(x, pars = "mu_beta_random", prob = prob)
        mu_beta_random <- .tidy_summary(mu_beta_random, c(gs$name, "predictor", "factor"), gs$map$label, pnames$location[P_random_ind], fnames)
        mu_beta_random <- .summary_rearrange(mu_beta_random, c(gs$name, "factor", "predictor"))
    }

    if(Q_random > 0) {
        logsd_beta_random <- .summarize(x, pars = "logsd_beta_random", prob = prob)
        logsd_beta_random <- .tidy_summary(logsd_beta_random, c(gs$name, "predictor", "factor"), gs$map$label, pnames$scale[Q_random_ind], fnames)
        logsd_beta_random <- .summary_rearrange(logsd_beta_random, c(gs$name, "factor", "predictor"))
    }

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
    PS <- x$meta$pred_spec
    pnames <- PS$pname
    fnames <- unlist(x$meta$indicator_spec$fname)

    mu_beta <- NA
    logsd_beta <- NA

    if(PS$P > 0) {
        mu_beta <- .summarize(x, pars = "mu_beta", prob = prob)
        mu_beta <- .tidy_summary(mu_beta, c("predictor", "factor"), pnames$location, fnames)
        mu_beta <- .summary_rearrange(mu_beta, c("factor", "predictor"))
    }

    if(PS$Q > 0) {
        logsd_beta <- .summarize(x, pars = "logsd_beta", prob = prob)
        logsd_beta <- .tidy_summary(logsd_beta, c("predictor", "factor"), pnames$scale, fnames)
        logsd_beta <- .summary_rearrange(logsd_beta, c("factor", "predictor"))
        
    }

    out <- nlist(mu_beta, logsd_beta)
    return(out)
}

.summary_between <- function(x, prob) {
    # Meta-data
    PS <- x$meta$pred_spec
    pnames <- PS$pname$between
    fnames <- unlist(x$meta$indicator_spec$fname)
    F <- x$meta$indicator_spec$F
    re_int_names <- paste0(fnames, "MAGICSEP", rep(c("mu", "logsd"), each = F))

    zeta <- NA
    if(PS$R > 0) {

        # Summarize
        zeta <- .summarize(x, pars = "zeta", prob = prob)

        # Tidy
        zeta <- .tidy_summary(zeta, c("predictor", "param"), pnames, re_int_names)

        # Separate
        zeta[, c("factor", "param")] <- .magicsep(zeta[, "param"], c("factor", "param"))

        # Rearrange
        zeta <- .summary_rearrange(zeta, c("factor", "param", "predictor"))
    }

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
##' @return x (Invisibly).
##' @author Stephen R. Martin
##' @export
print.summary.lmmelsm <- function(x, ...) {
    dots <- list(...)
    digits <- dots$digits %IfNull% x$meta$digits
    latent <- x$meta$latent
    facVarStr <- ifelse(latent, "Factor", "Variable") # To avoid some ifelses down the line

    # Diagnostics
    .sep()
    cat("Diagnostic Checks")
    .sep()
    .newline()
    .print.lmmelsm_diag(x$meta$stan$diag)
    .sep()

    IS <- x$meta$indicator_spec
    if(latent) {
        # Measurement model
        .sep()
        cat("Measurement Model")
        .sep()
        .newline()

        ## Loadings
        cat("Loadings")
        .newline()
        for(f in 1:IS$F) {
            cat("Factor: ")
            cat(IS$fname[[f]])
            .newline()
            with(x$summary, .print_table(lambda[lambda$factor == IS$fname[[f]] & lambda$item %in% IS$iname[[f]], ], digits, "factor"))
            .newline()
        }

        ## Residual SD
        .newline()
        cat("Residual Standard Deviation")
        .newline()
        with(x$summary, .print_table(sigma, digits))

        ## Intercepts
        .newline()
        cat("Intercepts")
        .newline()
        with(x$summary, .print_table(nu, digits))
    }
    
    # Factor Cors
    if(IS$F > 1) {
        .newline()
        cat(facVarStr, "correlations")
        .newline()
        with(x$summary, .print_table(factor_correlation, digits))
    }

    # Location model
    if(!latent) {
        .sep()
        cat("Location Intercepts")
        .sep()
        .newline()
        with(x$summary, .print_table(nu, digits))
    }
    if(x$meta$pred_spec$P > 0) {
        PS <- x$meta$pred_spec
        .sep()
        cat("Location model (Fixed effects)")
        .sep()
        .newline()
        for(f in 1:IS$F) {
            cat(facVarStr, ": ")
            cat(IS$fname[[f]])
            .newline()
            with(x$summary, .print_table(mu_coef[mu_coef$factor == IS$fname[[f]], ], digits, "factor"))
            .newline()
        }
    }

    # Scale model
    if(!latent) {
        .sep()
        cat("(Log) Scale Intercepts")
        .sep()
        .newline()
        with(x$summary, .print_table(sigma, digits))
    }
    if(x$meta$pred_spec$Q > 0) {
        PS <- x$meta$pred_spec
        .sep()
        cat("Scale model (Fixed effects)")
        .sep()
        .newline()
        for(f in 1:IS$F) {
            cat(facVarStr, ": ")
            cat(IS$fname[[f]])
            .newline()
            with(x$summary, .print_table(logsd_coef[logsd_coef$factor == IS$fname[[f]], ], digits, "factor"))
            .newline()
        }
    }

    # Between model
    if(x$meta$pred_spec$R > 0) {
        PS <- x$meta$pred_spec
        .sep()
        cat("Between-group scale model")
        .sep()
        .newline()
        for(f in 1:IS$F) {
            cat(facVarStr, ": ")
            cat(IS$fname[[f]])
            .newline()
            with(x$summary, .print_table(zeta[zeta$factor == IS$fname[[f]], ], digits, "factor"))
            .newline()
        }
    }

    # RE SDs
    .sep()
    cat("Random effect standard deviations")
    if(x$meta$pred_spec$R > 0) {
        .newline()
        .tab()
        cat("Note: Between-group scale model used.")
        .newline()
        .tab()
        cat("'mu' and 'logsd' represent RE-SDs when between-group covariates are zero.")
    }
    .sep()
    .newline()
    .print_table(x$summary$random_sigma, digits)

    # RE Cors
    ## Mean estimates only.
    .sep()
    cat("Random effect correlations (Posterior Means (Lower) and SDs (Upper))")
    .newline()
    .tab()
    cat("Note: See summary(out)$summary$Omega_mean_logsd for full summary.")
    .sep()
    .newline()

    cms <- x$summary$random_correlation

    re_total <- with(x$meta, indicator_spec$F*2 + pred_spec$P_random*indicator_spec$F + pred_spec$Q_random*indicator_spec$F)

    re_names <- .build_re_names(x$meta, sep = "_")
    
    corMat <- matrix(1, nrow = re_total, ncol = re_total)
    corMat[lower.tri(corMat)] <- x$summary$random_correlation$Mean
    corMat[upper.tri(corMat)] <- x$summary$random_correlation$SD
    rownames(corMat) <- colnames(corMat) <- re_names
    print(corMat, digits = digits)

    invisible(x)
}

.print_table <- function(x, digits, drop = NULL) {
    if(!is.null(drop)) {
        drop_cols <- match(drop, colnames(x))
        x <- x[, -drop_cols]
    }

    print.data.frame(x, digits = digits, row.names = FALSE)
    
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
##' @param prob Numeric (Default: .95; `[0 - 1]`). The desired mass to contain within the CrI.
##' @return Matrix.
##' @author Stephen R. Martin
##' @keywords internal
.summarize <- function(lmmelsm, pars, prob = .95) {
    if(inherits(lmmelsm, "lmmelsm")) {
        samps <- as.matrix(lmmelsm$fit, pars = pars)
    } else if (is.matrix(lmmelsm)) {
        samps <- lmmelsm
    }
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
##' @param x String. E.g., `"lambda[1,2]"`
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
##' .summarize creates an rstan-like summary with rownames, \code{mat[1:R, 1:C]}.
##' \code{.tidy_summary(mat, c("rows", "cols"))} would then create two new columns, "rows" and "cols" with the indices in them.
##' If arguments are provided in \code{...}, then these indicate the mappings between the indices and labeled values.
##' E.g., \code{.tidy_summary(mat, c("rows", "cols"), c("A", "B"), c("C", "D"))} would create two new columns, "rows" and "cols", and replace rows = 1 with rows = A; cols=2 with cols = D, and so on.
##' Useful for going from stan rownames, to labeled columns.
##' @title Takes stan summary, returns summary with indices-as-columns.
##' @param x Output of .summarize
##' @param labs The labels for each parameter index. E.g., "predictor", "factor"
##' @param ... Optional (but recommended). Mappings for indices. E.g., Index column 1 is replaced by ...`[[1]][col1Indices]`.
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
##' This function returns the stan-strings (when \code{string = TRUE}, e.g., \code{[2,1]}, \code{[3,1]}), or the row-index assuming column-major order.
##' @title Get indices for subsetting lower-tri summaries of square matrices. 
##' @param x Integer. Dimension of matrix.
##' @param string Logical (Default: FALSE). Whether to return strings (e.g., \code{[2,1]}, or row indices, assuming column-major ordering.)
##' @return Charactor vector (if \code{string} is TRUE) or integer vector.
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

.get_diagnostics <- function(object) {
    fit <- object$fit
    latent <- object$meta$latent

    bfmi <- rstan::get_bfmi(fit)
    bfmi_chains <- rstan::get_low_bfmi_chains(fit)

    div_num <- rstan::get_num_divergent(fit)
    div_iter <- rstan::get_divergent_iterations(fit)

    tree_num <- rstan::get_num_max_treedepth(fit)
    tree_iter <- rstan::get_max_treedepth_iterations(fit)

    # TODO: Limit this to parameters; not eta{_logsd}
    pars <- c("sigma", "nu",
              "mu_logsd_betas_random_sigma", "Omega_mean_logsd",
              "mu_random", "logsd_random", "mu_beta_random", "logsd_beta_random",
              "mu_beta", "logsd_beta",
              "zeta",
              "Omega_eta")
    if(latent) pars <- c(pars, "lambda")
    rhat <- rstan::summary(fit, pars = pars)$summary[,"Rhat"]
    rhat <- rhat[!is.na(rhat)]
    rhat_sorted <- sort(rhat, decreasing = TRUE)

    out <- nlist(bfmi,
                 bfmi_chains,
                 div_num,
                 div_iter,
                 tree_num,
                 tree_iter,
                 rhat,
                 rhat_sorted)
    class(out) <- "lmmelsm_diag"
    return(out)
}

.print.lmmelsm_diag <- function(x) {
    passed <- TRUE

    cat("Divergent Transitions: ")
    if(x$div_num > 0) {
        passed <- FALSE
        cat("Failed")
        .newline()
        .tab()
        cat(x$div_num, " divergences")
        .newline()
    } else {
        cat("Passed")
        .newline()
    }

    cat("Convergence: ")
    if(any(x$rhat >= 1.1)) {
        passed <- FALSE
        cat("Failed")
        .newline()
        .tab()
        cat("Ten highest Rhats: ")
        .newline()
        print(head(x$rhat_sorted, 10))
        .newline()
    } else {
        cat("Passed")
        .newline()
    }

    if(!passed) {
        cat("*** Diagnostics failed. Do not interpret estimates. ***")
        .newline()
    } else {
        cat("Diagnostics passed.")
        .newline()
    }

    # TODO: BFMI and treedepth?
}

.build_re_names <- function(meta, sep = "MAGICSEP") {
    pnames <- meta$pred_spec$pname
    fnames <- unlist(meta$indicator_spec$fname)
    F <- meta$indicator_spec$F
    P_random <- meta$pred_spec$P_random
    Q_random <- meta$pred_spec$Q_random
    P_random_ind <- meta$pred_spec$P_random_ind
    Q_random_ind <- meta$pred_spec$Q_random_ind

    re_names <- paste0(fnames, sep, rep(c("mu", "logsd"), each = F))

    if(P_random > 0) {
        re_names <- c(re_names, paste0(rep(fnames, each = P_random), sep, pnames$location[P_random_ind], "_mu"))
    }
    if(Q_random > 0) {
        re_names <- c(re_names, paste0(rep(fnames, each = Q_random), sep, pnames$scale[Q_random_ind], "_logsd"))
    }
    re_names
}
