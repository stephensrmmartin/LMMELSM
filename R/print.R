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
    out$summary$sigma <- .tidy_summary(out$summary$sigma, "item", ind_names)
    out$summary$nu <- .tidy_summary(out$summary$nu, "item", ind_names)

    # RE SDs.
    out$summary$mu_logsd_betas_random_sigma <- .summarize(object, pars = "mu_logsd_betas_random_sigma", prob = prob)
    out$summary$Omega_mean_logsd <- .summarize(object, pars = "Omega_mean_logsd", prob = prob)

    ## Restructure
    ### Get factor-param names for each RE.
    pnames <- out$meta$pred_spec$pname
    re_names <- paste0(fnames, "MAGICSEP", rep(c("mu", "logsd"), each = length(fnames)))
    re_total <- with(out$meta, 2 * indicator_spec$F + 2 * pred_spec$P_random + 2 * pred_spec$Q_random)
    if(out$meta$pred_spec$P_random > 0) { # If RE location slopes
        re_names <- c(re_names, paste0(rep(fnames, each = out$meta$pred_spec$P_random), "MAGICSEP", pnames$location[out$meta$pred_spec$P_random_ind]))
    }
    if(out$meta$pred_spec$Q_random > 0) { # IF RE scale slopes
        re_names <- c(re_names, paste0(rep(fnames, each = out$meta$pred_spec$Q_random), "MAGICSEP", pnames$scale[out$meta$pred_spec$Q_random_ind]))
    }
    re_int_names <- re_names[1:(2 * out$meta$indicator_spec$F)]
    re_slope_names <- re_names[(2 * out$meta$indicator_spec$F + 1):re_total]

    ### Restructure them
    out$summary$mu_logsd_betas_random_sigma <- .tidy_summary(out$summary$mu_logsd_betas_random_sigma, "param", re_names)
    out$summary$Omega_mean_logsd <- .tidy_summary(out$summary$Omega_mean_logsd, c("row", "col"), re_names, re_names)
    ## Split MAGICSEP'd names into columns.
    out$summary$mu_logsd_betas_random_sigma[, c("factor", "param")] <- .magicsep(out$summary$mu_logsd_betas_random_sigma[, "param"], c("factor", "param"))

    out$summary$Omega_mean_logsd[, c("row_factor", "row_param")] <- .magicsep(out$summary$Omega_mean_logsd[, "row"], c("row_factor", "row_param"))
    out$summary$Omega_mean_logsd[, c("col_factor", "col_param")] <- .magicsep(out$summary$Omega_mean_logsd[, "col"], c("col_factor", "col_param"))
    out$summary$Omega_mean_logsd[, c("row", "col")] <- NULL
    out$summary$Omega_mean_logsd <- out$summary$Omega_mean_logsd[.full_to_lower_tri(re_total), ]

    # Location predictors
    out$summary$mu_beta <- .summarize(object, pars = "mu_beta", prob = prob)
    out$summary$mu_beta <- .tidy_summary(out$summary$mu_beta, c("predictor", "factor"), pnames$location, fnames)

    # Scale predictors
    out$summary$logsd_beta <- .summarize(object, pars = "logsd_beta", prob = prob)
    out$summary$logsd_beta <- .tidy_summary(out$summary$logsd_beta, c("predictor", "factor"), pnames$scale, fnames)

    # L2 scale predictors
    out$summary$zeta <- .summarize(object, pars = "zeta", prob = prob)
    out$summary$zeta <- .tidy_summary(out$summary$zeta, c("predictor", "param"), pnames$between, re_int_names) # RE intercepts only
    out$summary$zeta[, c("factor", "param")] <- .magicsep(out$summary$zeta[, "param"], c("factor", "param"))

    # Ranefs
    out$summary$mu_random <- .summarize(object, pars = "mu_random", prob = prob)
    out$summary$logsd_random <- .summarize(object, pars = "logsd_random", prob = prob)
    out$summary$mu_beta_random <- .summarize(object, pars = "mu_beta_random", prob = prob)
    out$summary$logsd_beta_random <- .summarize(object, pars = "logsd_beta_random", prob = prob)

    out$summary$mu_random <- .tidy_summary(out$summary$mu_random, c(out$meta$group_spec$name, "factor"), out$meta$group_spec$map$label, fnames)
    out$summary$logsd_random <- .tidy_summary(out$summary$logsd_random, c(out$meta$group_spec$name, "factor"), out$meta$group_spec$map$label, fnames)
    out$summary$mu_beta_random <- .tidy_summary(out$summary$mu_beta_random, c(out$meta$group_spec$name, "predictor", "factor"), out$meta$group_spec$map$label, pnames$location[out$meta$pred_spec$P_random_ind], fnames)
    out$summary$logsd_beta_random <- .tidy_summary(out$summary$logsd_beta_random, c(out$meta$group_spec$name, "predictor", "factor"), out$meta$group_spec$map$label, pnames$scale[out$meta$pred_spec$Q_random_ind], fnames)

    # Eta correlations
    out$summary$Omega_eta <- .summarize(object, pars = "Omega_eta", prob = prob)
    out$summary$Omega_eta <- .tidy_summary(out$summary$Omega_eta, c("row", "col"), fnames, fnames)
    out$summary$Omega_eta <- out$summary$Omega_eta[.full_to_lower_tri(out$meta$indicator_spec$F), ]

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
