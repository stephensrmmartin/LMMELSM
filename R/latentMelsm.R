##' Fits a MELSM on latent variables.
##'
##' Currently supports multiple latent factors, exogenous, and exogenous variables.
##' Data are assumed to be two-level data. I.e., multiple indicators, repeatedly measured within person.
##' Currently assumes measurement invariance (i.e., the measurement model params are equivalent across groups).
##' Excludes rows with missing data.
##' The Stan model currently uses the unit-variance identification.
##' @title Fit two-level latent MELSM.
##' @param formula Formula or list of formulas. LHS of each should be factor name, RHS should be indicators.
##' @param group Raw grouping variable name (not character).
##' @param data Data frame.
##' @param ... Options passed onto \code{\link[rstan]{sampling}}
##' @return melsm_latent object.
##' @author Stephen R. Martin
##' @import rstan
melsm_latent <- function(formula, group, data, ...) {
    # TODO: Add default sampling handling
    # TODO: For a package, make sure to use stanmodels list instead.

    d <- .parse_formula(formula, group = substitute(group), data)

    sm <- stan_model("../Stan/melsm_2l_mi.stan")
    ## sm <- stan_model("../Stan/melsm_2l_mi_glm.stan")
    pars <- c("nu",
              "lambda",
              "sigma",
              "eta_mean",
              "eta_sd",
              "eta",
              "sigma_mean_logsd",
              "Omega_eta",
              "Omega_mean_logsd")
    sOut <- sampling(sm,
                    data = d$stan_data,
                    cores = 4,
                    iter = 1000,
                    pars = pars,
                    control = list(adapt_delta = .95)
                    )
    return(sOut)

}
##' @title Convert spec to stan data.
##' @param formulaList Formula or list of formulas.
##' @param group Group symbol.
##' @param data Data frame.
##' @return List.
##' @author Stephen R. Martin
##' @import Formula
.parse_formula <- function(formulaList, group, data) {
    # TODO: Allow exogenous predictors; endogenous outcomes.
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
    lhs_names <- lapply(flist, function(f){
       all.vars(f)[1]
    })
    rhs_names <- lapply(flist, function(f) {
        if(formula) {
            attr(terms(f), "term.labels")
        } else {
            all.vars(f)[-1]
        }
    })
    out <- list(factor = lhs_names, indicator = rhs_names)

    return(out)
}

##' @title Combines multiple formulas' RHS into one.
##' @param flist 
##' @return RHS formula.
##' @author Stephen R. Martin
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
