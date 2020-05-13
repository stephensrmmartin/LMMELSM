##' Fits a MELSM on latent variables.
##'
##' Currently supports multiple latent factors, exogenous, and exogenous variables.
##' Data are assumed to be two-level data. I.e., multiple indicators, repeatedly measured within person.
##' Currently assumes measurement invariance (i.e., the measurement model params are equivalent across groups).
##' Excludes rows with missing data.
##' The Stan model currently uses the unit-variance identification.
##'
##' Formulas should be of the form:
##'
##' \code{latentName ~ indicator1 + indicator2 + indicator3},
##' where \code{latentName} is the name of the latent variable.
##' The right-hand side (RHS) should contain the variables in \code{data} onto which the factor loads.
##' In the simplest case, only one formula is required. This defines a single-factor model.
##' If multiple latent factors are desired, a list of formulas can be provided.
##'
##' To predict the latent locations and scales, include the respective formulas in the list of formulas.
##' For example:
##' 
##' \code{location ~ x1 + x2}
##' \code{scale ~ x1 + x2}
##'
##' The covariates do not need to be the same across the location and scale.
##' The specified covariates will be used to predict the location and scale of \emph{all} latent factors via multivariate regression.
##'
##' \emph{Note}: Because \code{location} and \code{scale} represent special formulas, latent factors cannot be called location or scale.
##' It is assumed that any formula with \code{location} or \code{scale} on the left-hand side (LHS) is a predictive formula, not a latent variable specification.
##' @title Fit two-level latent MELSM.
##' @param formula Formula or list of formulas. LHS of each should be factor name, RHS should be indicators.
##' @param group Raw grouping variable name (not character).
##' @param data Data frame.
##' @param ... Options passed onto \code{\link[rstan]{sampling}}
##' @return melsm_latent object.
##' @author Stephen R. Martin
##' @import rstan
##' @importFrom parallel detectCores
melsm_latent <- function(formula, group, data, ...) {
    # TODO: Add default sampling handling
    # TODO: For a package, make sure to use stanmodels list instead.
    # Set defaults
    dots <- list(...)
    stan_args <- list()
    stan_args$control <- dots$control %IfNull% list(adapt_delta = .95)
    stan_args$control$adapt_delta <- stan_args$control$adapt_delta %IfNull% .95
    stan_args$cores <- dots$cores %IfNull% detectCores()
    stan_args$chains <- dots$chains %IfNull% 4
    stan_args$iter <- dots$iter %IfNull% 2000
    ## Remove from dots the things that are specified here
    dots[names(dots) %in% names(stan_args)] <- NULL

    d <- .parse_formula(formula, group = substitute(group), data)

    stan_args$model <- stanmodels$melsm2MIGLM

    pars <- c("nu",
              "lambda",
              "sigma",
              "eta_mean",
              "eta_sd",
              "eta",
              "sigma_mean_logsd",
              "Omega_eta",
              "Omega_mean_logsd")
    stan_args$pars <- pars

    sOut <- do.call(sampling, c(stan_args, dots))

    return(sOut)

}

##' @title Convert spec to stan data.
##' @param formula Formula. Should be either location ~ ... or scale ~ ... .
##' @param group Raw group name.
##' @param data Data frame.
##' @return List.
##' @author Stephen R. Martin
.parse_prediction_formula <- function(formula, group, data) {
    # TODO: Implement this
    
}

##' @title Convert spec to stan data.
##' @param formulaList Formula or list of formulas.
##' @param group Group symbol.
##' @param data Data frame.
##' @return List.
##' @author Stephen R. Martin
##' @import Formula
##' @keywords internal
.parse_formula2 <- function(formulaList, group, data) {
    # Make it a list of formulas
    if(!is.list(formulaList)) {
        flist <- list(formulaList)
    } else {
        flist <- formulaList
    }

    # Output structure
    out <- list(meta = list(), stan_data = list())

    # Convert to Formula::Formula
    flist <- lapply(flist, as.Formula)

    # Separate location/scale formulas from measurement formulas.
    ## mlist: List of measurement fact formulas
    ## plist: List of predictive location/scale formulas
    ## flist: All formulas (good for model framing)
    which_loc_sca <- .which_location_scale(flist, reduce = TRUE)
    plist <- flist[which_loc_sca]
    names(plist) <- sapply(plist, .get_LHS)
    mlist <- flist
    mlist[which_loc_sca] <- NULL

    # Check for LHS names
    for(f in seq_len(length(flist))) {
        if(length(flist[[f]])[1] != 1) {
            stop("Factor name should be provided on LHS of formula.")
        }
    }

    # Get group spec
    group <- substitute(group)
    group_name <- deparse(group)

    # Create model frame (All predictors, indicators, and the grouping variable)
    flist_RHS <- .combine_RHS(flist)
    mf <- model.frame(flist_RHS, data, na.action = na.pass)
    mf[,group_name] = data[,group_name]

    # Remove missings
    removed_ind <- which(!complete.cases(mf)) 
    removed_N <- length(removed_ind)
    out$meta$missings = list(N = removed_N, ind = removed_ind)

    mf <- mf[-removed_ind,]
    if(removed_N >= 1) {
        warning("Removing", removed_N, "incomplete cases.")
    }
    out$data = mf

    # Group object
    group <- list(name = group_name,
                  data = mf[, group_name],
                  numeric = as.numeric(as.factor(mf[, group_name])))
    out$meta$group <- group
    out$stan_data$group <- group$numeric


    # Get indicator matrix
    mlist_RHS <- .combine_RHS(mlist)
    mm <- model.matrix(mlist_RHS, mf)[, -1] # No intercept
    ind_spec <- .get_indicator_spec(mm, mlist)
    out$stan_data[c("J_f", "F_ind")] <- ind_spec[c("J_f", "F_ind")]
    out$stan_data$y <- mm

    # Predictor matrices
    pred_spec <- .parse_formula.predictor(plist, mf, group$data)
}

.parse_formula.predictor <- function(plist, mf, group) {
    plist$location <- plist$location %IfNull% Formula(location ~ 1)
    plist$scale <- plist$scale %IfNull% Formula(scale ~ 1)
    mf <- model.frame(.combine_RHS(plist), mf)

    x_loc <- model.matrix(plist$location, mf)[,-1]
    x_sca <- model.matrix(plist$scale, mf)[,-1]
    P <- ncol(x_loc)
    Q <- ncol(x_sca)

    if(length(plist$location)[2] == 2) {
        P_random_RHS <- .get_RHS(formula(plist$location, rhs = 2))
        P_random <- length(P_random_RHS)
        P_random_ind <- match(P_random_RHS, colnames(x_loc))
    } else {
        P_random <- 0
        P_random_ind <- integer()
    }

    if(length(plist$scale)[2] == 2) {
        Q_random_RHS <- .get_RHS(formula(plist$scale, rhs = 2))
        Q_random <- length(Q_random_RHS)
        Q_random_ind <- match(Q_random_RHS, colnames(x_sca))
    } else {
        Q_random <- 0
        Q_random_ind <- integer()
    }

    # Specify efficiency options
    intercept_only <- P == 0 & Q == 0
    L2_pred_only <- .detect_L2_only(mf, group)


    out <- nlist(intercept_only,
                 L2_pred_only,
                 P,
                 Q,
                 P_random,
                 Q_random,
                 P_random_ind,
                 Q_random_ind,
                 x_loc,
                 x_sca
                 )
    return(out)
}

##' @title Convert spec to stan data.
##' @param formulaList Formula or list of formulas.
##' @param group Group symbol.
##' @param data Data frame.
##' @return List.
##' @author Stephen R. Martin
##' @import Formula
##' @keywords internal
.parse_formula <- function(formulaList, group, data) {
    # TODO: Allow exogenous predictors; endogenous outcomes.
    # TODO: Allow 2-level predictive formulas to separate out L1/L2 (random/fixed), or just random/fixed.
    # TODO: Redo this function from the ground-up
    if(!is.list(formulaList)) {
        forms <- list(formulaList)
    } else {
        forms <- formulaList
    }

    forms <- lapply(forms, as.Formula)

    # Check for LHS
    for(f in seq_len(length(forms))) {
        if(length(forms[[f]])[1] != 1) {
            stop("Factor name should be provided on LHS of formula.")
        }
    }

    # Predictor formulas
    ## Which formulas correspond to location and scale
    which_loc_sca <- .which_location_scale(forms)

    # Get LHS and RHS terms
    fterms <- .get_formula_names(forms, formula = TRUE)

    # Make model.frame
    RHS.inds <- .combine_RHS(forms)
    mf.inds <- model.frame(RHS.inds, data = data, na.action = na.pass)

    # Grouping data
    group_name <- deparse(group)
    group_orig <- data[[group_name]]


    # Remove missings
    which_na <- lapply(cbind(mf.inds, data[[group_name]]), function(x){
        which(is.na(x))
    })
    which_na_vec <- unique(do.call(c, which_na))

    if(length(which_na_vec) >= 1) {
        warning("Removing ", length(which_na_vec), "rows with NA values.")
        mf.inds <- mf.inds[-which_na_vec,]
        group_orig <- group_orig[-which_na_vec]
    }

    group_numeric <- as.numeric(as.factor(group_orig))
    group_K <- length(unique(group_numeric))
    group_data <- list(name = group_name, group = group_orig, group_code = group_numeric, group_K = group_K)

    # Make model.matrix
    mm.inds <- model.matrix(RHS.inds, data = mf.inds)[, -1] ## No intercept

    # Get indicator spec
    ind.spec <- .get_indicator_spec(mm.inds, forms)

    # Meta-data
    N <- nrow(mf.inds)
    J <- ncol(mm.inds)
    `F` <- length(fterms$factor)

    # Package up
    mf <- mf.inds
    mf[[group_name]] <- group_orig

    stan_data <- list(N = N,
                      J = J,
                      `F` = `F`,
                      K = group_K,
                      group = group_numeric,
                      J_f = ind.spec$J_f,
                      F_ind = ind.spec$F_ind,
                      x = mm.inds)
    out <- list(meta = list(N = N,
                            J = J,
                            `F` = `F`,
                            K = group_K,
                            group_data = group_data,
                            formula = forms,
                            fnames = fterms$factor,
                            inames = fterms$indicator),
                data = mf,
                stan_data = stan_data)

    return(out)

}

##' @title Get names in formula.
##' @param flist List of formulas.
##' @param formula Logical. Whether to return the raw RHS (TRUE) or the vars needed (FALSE).
##' @return List of lists.
##' @author Stephen R. Martin
.get_formula_names <- function(flist, formula = TRUE) {
    lhs_names <- lapply(flist, .get_LHS)
    rhs_names <- lapply(flist, .get_RHS, terms = formula)

    out <- list(factor = lhs_names, indicator = rhs_names)

    return(out)
}

##' @title Get LHS variable as string.
##' @param formula Formula.
##' @return String.
##' @author Stephen R. Martin
##' @keywords internal
.get_LHS <- function(formula) {
    lhs_name <- all.vars(formula)[1]

    return(lhs_name)
}
##' @title Get RHS terms or variables.
##' @param formula Formula.
##' @param terms Logical (Default: TRUE). Whether to return the terms (TRUE) or the variables (FALSE). E.g., "I(x^2)" vs "x"
##' @return Character vector.
##' @author Stephen R. Martin
##' @keywords internal
.get_RHS <- function(formula, terms = TRUE) {
    rhs_name <- ifelse(terms, attr(terms(formula), "term.labels"), all.vars(formula)[-1])

    return(rhs_name)
}

##' @title Combines multiple formulas' RHS into one.
##' @param flist 
##' @return RHS formula.
##' @author Stephen R. Martin
##' @keywords internal
.combine_RHS <- function(flist) {
    fnames <- .get_formula_names(flist, formula = TRUE)
    RHS <- do.call(c, fnames$indicator)
    RHS_combined <- as.Formula(paste0("~ ", paste0(unique(RHS), collapse = " + ")))
    return(RHS_combined)
}

##' @title Get indicator spec for stan model.
##' @param mm.inds Model matrix. Complete data of all non-duplicated indicators for all factors.
##' @param flist Formula list.
##' @return List with J_f (vector; Number of Indicators for each factor) and F_ind (matrix; each row contains the columns in mm.inds for the row-th factor).
##' @author Stephen R. Martin
##' @keywords internal
.get_indicator_spec <- function(mm.inds, flist) {
    fnames <- .get_formula_names(flist, formula = TRUE)
    mmnames <- colnames(mm.inds)
    `F` <- length(flist)
    J <- ncol(mm.inds)

    J_f <- array(0, dim = `F`)
    F_ind <- matrix(0, `F`, J)

    for(f in seq_len(`F`)) {
        J_f[f] <- length(fnames[["indicator"]][[f]])
        F_ind[f, seq_len(J_f[f])] <- match(fnames[["indicator"]][[f]], mmnames)
    }

    out <- list(J_f = J_f, F_ind = F_ind)
    return(out)
}

##' @title Check for location-scale formulas
##' @param flist Formula list.
##' @return Numeric vector (length 2); which formulas in flist correspond to location, scale.
##' @author Stephen R. Martin
##' @keywords internal
.which_location_scale <- function(flist, reduce = TRUE) {
    lhs_names <- tolower(lapply(flist, .get_LHS))
    loc_scale <- match(c("location", "scale"), lhs_names)
    names(loc_scale) <- c("location", "scale")

    # Check whether multiple location/scales exist
    which_lhs_match_location <- lhs_names %in% c("location")
    which_lhs_match_scale <- lhs_names %in% c("scale")
    if(sum(which_lhs_match_location) > 1) {
        stop("Multiple formulas for 'location' provided.")
    }
    if(sum(which_lhs_match_scale) > 1) {
        stop("Multiple formulas for 'scale' provided.")
    }

    if(reduce) {
        loc_scale <- loc_scale[!is.na(loc_scale)]
    }

    return(loc_scale)
}

##' Detects whether the predictors are level-2 only.
##' E.g., whether location and scale are covariates are the same across all n_k observations of each K.
##' This is important for efficiency reasons.
##' If the covariates are invariant across repeated observations of the given group k, for all K, then we can compute predicted values once, and broadcast the prediction, rather than compute the prediction for every single row.
##' Specifically, it detects if all x == x[1], where x is a group's data, for each column in mf.
##' @title Detect whether the predictors are L2-only
##' @param mf Data frame for predictors. Should contain no missings.
##' @param group Grouping variable for the model frame.
##' @return Logical. TRUE if the covariates appear to be level-2 only.
##' @author Stephen R. Martin
##' @keywords internal
.detect_L2_only <- function(mf, group) {
    col_same <- sapply(mf, function(col) {
        group_same <- tapply(col, group, function(x) {
            all(x == x[1])
        })
        all(group_same)
    })
    return(all(col_same))
}
