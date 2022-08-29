##' Fits a mixed effects location scale model on one or more observed or latent variables.
##' Currently supports multiple endogenous latent factors or observed outcomes,  and exogenous observed variables.
##' Data are assumed to be two-level data. I.e., multiple indicators, repeatedly measured within group.
##' Currently assumes measurement invariance (i.e., the measurement model params are equivalent across groups) and a unit-variance identification for latent variables.
##' Excludes rows with missing data (and warns the user).
##'
##' @section Model specification:
##' 
##' The model is specified as a list of formulas.
##' LMMELSM supports the specification of latent measurement models, location models, scale models, between-group scale models, and (if latent variables are undesired) observed outcome models.
##' The covariates do not need to be the same across the location, scale, and between-group models.
##' The specified covariates will be used to predict the location and scale of \emph{all} latent factors via multivariate regression.
##' 
##' The latent factor model is specified as follows.
##' In the simplest case, only one formula is required, and a single-factor model is estimated.
##' The left-hand side (LHS) contains the user-specified latent variable name, and the right-hand side (RHS) contains the indicators.
##' Let "latent1" and "latent2" be user-chosen names of two latent variables with three indicators each.
##' Then the formula syntax would be:
##' \code{list(latent1 ~ y1 + y2 + y3, latent2 ~ y4 + y5 + y6)}
##'
##'
##' The location model is specified as either a one or two-part formula.
##' The LHS must be "location" and the RHS contains the covariates.
##' Random slopes are specified in the optional second part, separated by "|".
##' Because LMMELSM fits MELSMs, random intercepts are \emph{always} included.
##' For example, if x1 and x2 are two location predictors, then:
##' 
##' \code{location ~ x1 + x2}
##'
##' specifies a location model with a random intercept per factor, and no random slopes.
##'
##' \code{location ~ x1 + x2 | x1}
##'
##' specifies a location model with a random intercept per factor, a random x1 coefficient per factor, and no random x2 coefficient.
##'
##' The within-group scale model is specified similarly.
##' The LHS must be "scale" and the RHS contains the covariates.
##' Random intercepts are always included, and random slopes are specified in the optional second part of the RHS.
##' For example, if x2 and x3 are two scale predictors, then:
##'
##' \code{scale ~ x2 + x3}
##'
##' specifies a scale model with a random intercept per factor, and no random slopes.
##'
##' \code{scale ~ x2 + x3 | x3}
##'
##' specifies a scale model with a random intercept perfactor, a random x3 coefficient per factor, and no random x2 coefficient.
##'
##' The between-group scale model is specified by a LHS of "between" and RHS containing covariates.
##' There are no random coefficients permitted in the between-group scale model.
##' The between-group scale model is responsible for modeling the random effect standard deviations.
##' \emph{Note:} The between-group model \emph{only} models the SDs of the random location and scale \emph{intercepts}.
##'
##' \code{between ~ x2}
##'
##' specifies a between-group scale model on the SDs of the location and scale intercepts for each factor.
##'
##' If you want to fit a non-latent multivariate MELSM, use "observed" as the LHS:
##'
##' For example, if y1, y2, and y3 are three observed outcome variables, then
##' 
##' \code{observed ~ y1 + y2 + y3}
##'
##' would fit an M-MELSM.
##' Location, scale, and between-group models can still be specified, but they will model the observed variables, rather than latent variables.
##' You cannot currently have both observed and latent outcomes in the same model.
##'
##' \emph{Note}: Because \code{location}, \code{scale}, \code{between}, and \code{observed} represent special formulas, latent factors cannot be named location, scale, between, nor observed.
##' It is assumed that any formula with \code{location}, \code{scale}, or \code{between} on the left-hand side (LHS) is a predictive formula, not a latent variable specification.
##' @title Specify and fit the (latent) (multivariate) melsm.
##' @param formula Formula or list of formulas. See section on model specification.
##' @param group Raw grouping variable name (not character).
##' @param data Data frame.
##' @param ... Options passed onto \code{\link[rstan]{sampling}}
##' @return lmmelsm object.
##' @author Stephen R. Martin
##' @rawNamespace import(rstan, except = loo)
##' @importFrom parallel detectCores
##' @importFrom stats complete.cases dnorm formula model.frame model.matrix na.pass quantile rnorm
##' @importFrom utils head strcapture
##' @export
##' @examples
##' \donttest{
##' data(sim_data)
##'
##' # Fit LMMELSM with two latent factors (A and B),
##' # Location model with one random coefficient
##' # Scale model with one random coefficient
##' # Between-group scale model with one covariate
##' fit <- lmmelsm(list(A ~ A_1 + A_2 + A_3 + A_4 + A_5 + A_6,
##'                     B ~ N_1 + N_2 + N_3 + N_4 + N_5 + N_6,
##'                     location ~ x1 + baseline | x1,
##'                     scale ~ x2 + baseline | x2,
##'                     between ~ baseline),
##'                subject, sim_data, cores = 2, iter = 500, chains = 2
##'               )
##'
##' # Summarize fit
##' summary(fit)
##'
##' # Get random effects
##' ranef(fit)
##' # Get group-specific parameter values
##' coef(fit)
##' # Get approximate leave-one-out
##' loo(fit)
##' }
lmmelsm <- function(formula, group, data, ...) {
    # Set defaults
    dots <- list(...)
    stan_args <- list()
    stan_args$control <- dots$control %IfNull% list(adapt_delta = .95)
    stan_args$control$adapt_delta <- stan_args$control$adapt_delta %IfNull% .95
    stan_args$cores <- dots$cores %IfNull% detectCores()
    stan_args$chains <- dots$chains %IfNull% 4
    stan_args$iter <- dots$iter %IfNull% 2000
    stan_args$init <- dots$init %IfNull% 0
    ## Remove from dots the things that are specified here
    dots[names(dots) %in% names(stan_args)] <- NULL

    d <- .parse_formula(formula, group = substitute(group), data)

    stan_args$object <- stanmodels$lmmelsmPred
    if(!d$meta$latent) {
        stan_args$object <- stanmodels$lmmelsmPredObs2
    }
    stan_args$data <- d$stan_data
    stan_args$data$prior_only <- dots$prior_only %IfNull% FALSE
    dots$prior_only <- NULL

    pars <- c("nu",
              "lambda",
              "sigma",
              "eta",
              "eta_logsd",
              "mu_beta",
              "logsd_beta",
              "zeta",
              "mu_random",
              "logsd_random",
              "mu_beta_random",
              "logsd_beta_random",
              "mu_logsd_betas_random_sigma",
              "Omega_eta",
              "Omega_mean_logsd")
    if(!d$meta$latent) {
        pars <- pars[-2] # Remove lambda
    }
    stan_args$pars <- pars

    sOut <- suppressWarnings(do.call(sampling, c(stan_args, dots)))

    out <- list(fit = sOut,
                meta = d$meta,
                data = d$data,
                stan_data = d$stan_data,
                stan_args = stan_args
                )

    class(out) <- "lmmelsm"

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
    names(flist) <- sapply(flist, .get_LHS)

    # Detected observedness
    observed <- FALSE
    if("observed" %in% tolower(names(flist))) {
        observed <- TRUE
    }

    # Separate location/scale formulas from measurement formulas.
    ## mlist: List of measurement factor formulas
    ## plist: List of predictive location/scale/between formulas
    ## flist: All formulas (good for model framing)
    plist <- list(location = Formula(location ~ 1), scale = Formula(scale ~ 1), between = Formula(between ~ 1))
    which_loc_sca <- .which_location_scale(flist, reduce = TRUE)
    plist[names(which_loc_sca)] <- flist[which_loc_sca]
    mlist <- flist
    mlist[which_loc_sca] <- NULL

    # Check for LHS names
    for(f in seq_len(length(flist))) {
        if(length(flist[[f]])[1] != 1) {
            stop("Factor name should be provided on LHS of formula.")
        }
    }

    # Get group spec
    group <- group
    group_name <- deparse(group)

    # Create model frame (All predictors, indicators, and the grouping variable)
    flist_RHS <- .combine_RHS(flist)
    mf <- model.frame(flist_RHS, data, na.action = na.pass)
    mf[,group_name] = data[,group_name]

    # Remove missings
    removed_ind <- which(!complete.cases(mf)) 
    removed_N <- length(removed_ind)
    out$meta$missings = list(N = removed_N, ind = removed_ind)

    if(removed_N >= 1) {
        mf <- mf[-removed_ind,]
        warning("Removing ", removed_N, " incomplete cases.")
    }
    out$data = mf

    # Group object
    group_spec <- list(name = group_name,
                       data = mf[, group_name],
                       numeric = as.numeric(as.factor(mf[, group_name])),
                       K = length(unique(mf[, group_name])))
    group_spec$map <- data.frame(numeric = 1:group_spec$K,
                                 label = group_spec$data[match(1:group_spec$K, group_spec$numeric)])
    out$meta$group_spec <- group_spec
    out$stan_data$group <- group_spec$numeric
    out$stan_data$K <- group_spec$K


    # Get indicator matrix
    if(observed) {
        indicator_spec <- .parse_formula.observed(mlist, mf)
        mlistNames <- .get_formula_names(mlist, formula = TRUE)
        out$meta$indicator_spec <- indicator_spec
        out$meta$latent <- FALSE
        out$meta$indicator_spec$fname <- as.list(mlistNames$indicator$observed)
        out$meta$indicator_spec$iname <- mlistNames$indicator
        out$meta$indicator_spec$mname <- mlistNames$indicator$observed
        out$meta$indicator_spec$mlist <- mlist
        out$stan_data <- c(out$stan_data, indicator_spec)
    } else {
        indicator_spec <- .parse_formula.indicators(mlist, mf)
        mlistNames <- .get_formula_names(mlist, formula = TRUE)
        out$meta$latent <- TRUE
        out$meta$indicator_spec <- indicator_spec
        out$meta$indicator_spec$fname <- mlistNames$factor
        out$meta$indicator_spec$iname <- mlistNames$indicator
        out$meta$indicator_spec$mname <- colnames(indicator_spec$y)
        out$meta$indicator_spec$mlist <- mlist
        out$stan_data <- c(out$stan_data, indicator_spec)
    }

    # Predictor matrices
    pred_spec <- .parse_formula.predictor(plist, mf, group_spec$data)
    plistNames <- .get_formula_names(plist)
    out$meta$pred_spec <- pred_spec
    out$meta$pred_spec$pname <- plistNames$indicator
    out$meta$pred_spec$plist <- plist
    out$stan_data <- c(out$stan_data, pred_spec)

    # Misc
    out$meta$formula <- flist

    return(out)
}

##' @title Compute indicator data.
##' @param mlist List of measurement formulas.
##' @param mf Data frame with complete cases.
##' @return Named List.
##' @author Stephen R. Martin
##' @keywords internal
.parse_formula.indicators <- function(mlist, mf) {
    mlist_RHS <- .combine_RHS(mlist)
    mm <- model.matrix(mlist_RHS, mf)[, -1] # No intercept
    ind_spec <- .get_indicator_spec(mm, mlist)

    out <- nlist(y = mm,
                 J_f = ind_spec$J_f,
                 F_ind = ind_spec$F_ind,
                 J = ncol(mm),
                 F = length(mlist),
                 N = nrow(mm)
                 )
    return(out)
}

.parse_formula.observed <- function(mlist, mf) {
    mlist_RHS <- .combine_RHS(mlist)
    mm <- model.matrix(mlist_RHS, mf)[, -1, drop = FALSE]

    J <- ncol(mm)
    `F` <- J
    J_f <- as.array(rep(1, `F`))
    F_ind <- matrix(0, `F`, J)
    F_ind[,1] <- seq_len(J)

    out <- nlist(y = mm,
                 J_f,
                 F_ind,
                 J,
                 `F`,
                 N = nrow(mm))

    return(out)
}

##' @title Compute predictor data.
##' @param plist List of prediction formulas.
##' @param mf Data frame containing complete cases.
##' @param group Grouping data.
##' @return Named list.
##' @author Stephen R. Martin
##' @keywords internal
.parse_formula.predictor <- function(plist, mf, group) {
    ## plist$location <- plist$location %IfNull% Formula(location ~ 1)
    ## plist$scale <- plist$scale %IfNull% Formula(scale ~ 1)
    mf <- model.frame(.combine_RHS(plist), mf)

    x_loc <- model.matrix(plist$location, mf)[,-1, drop = FALSE]
    x_sca <- model.matrix(plist$scale, mf)[,-1, drop = FALSE]
    x_bet <- model.matrix(plist$between, mf)[, -1, drop = FALSE]
    P <- ncol(x_loc)
    Q <- ncol(x_sca)
    R <- ncol(x_bet)

    if(length(plist$location)[2] == 2) {
        P_random_RHS <- .get_RHS(formula(plist$location, rhs = 2))
        P_random <- length(P_random_RHS)
        P_random_ind <- match(P_random_RHS, colnames(x_loc))
        P_random_ind <- array(P_random_ind, dim = c(P_random))
    } else {
        P_random <- 0
        P_random_ind <- integer()
    }

    if(length(plist$scale)[2] == 2) {
        Q_random_RHS <- .get_RHS(formula(plist$scale, rhs = 2))
        Q_random <- length(Q_random_RHS)
        Q_random_ind <- match(Q_random_RHS, colnames(x_sca))
        Q_random_ind <- array(Q_random_ind, dim = c(Q_random))
    } else {
        Q_random <- 0
        Q_random_ind <- integer()
    }

    # Specify efficiency options
    ## intercept_only <- P == 0 & Q == 0
    mf_loc_sca <- model.frame(.combine_RHS(plist[c("location","scale")]), mf)
    L2_pred_only <- .detect_L2_only(mf_loc_sca, group)


    out <- nlist(L2_pred_only,
                 P,
                 Q,
                 R,
                 P_random,
                 Q_random,
                 P_random_ind,
                 Q_random_ind,
                 x_loc,
                 x_sca,
                 x_bet
                 )
    return(out)
}

##' @title Get names in formula.
##' @param flist List of formulas.
##' @param formula Logical. Whether to return the raw RHS (TRUE) or the vars needed (FALSE).
##' @return List of lists.
##' @author Stephen R. Martin
##' @keywords internal
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
    if(terms) {
        out <- attr(terms(formula), "term.labels")
        if(length(out) == 0) {
            out <- "1" # Edge case: If intercept-only, then return 1.
        }
        return(out)
    } else {
        out <- all.vars(formula)[-1]
        if(length(out) == 0) {
            out <- "1"
        }
        return(out)
    }
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
    loc_scale <- match(c("location", "scale", "between"), lhs_names)
    names(loc_scale) <- c("location", "scale", "between")

    # Check whether multiple location/scales exist
    which_lhs_match_location <- lhs_names %in% c("location")
    which_lhs_match_scale <- lhs_names %in% c("scale")
    which_lhs_match_between <- lhs_names %in% c("between")
    if(sum(which_lhs_match_location) > 1) {
        stop("Multiple formulas for 'location' provided.")
    }
    if(sum(which_lhs_match_scale) > 1) {
        stop("Multiple formulas for 'scale' provided.")
    }
    if(sum(which_lhs_match_between) > 1) {
        stop("Multiple formulas for 'between' provided.")
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
##' Specifically, it detects if all \code{x == x[1]}, where x is a group's data, for each column in mf.
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

has_latent <- function(lmmelsm) {
    lmmelsm$meta$latent
}

has_location <- function(lmmelsm) {
    lmmelsm$meta$pred_spec$P > 0
}

has_scale <- function(lmmelsm) {
    lmmelsm$meta$pred_spec$Q > 0
}

has_between <- function(lmmelsm) {
    lmmelsm$meta$pred_spec$R > 0
}

has_random_location <- function(lmmelsm) {
    lmmelsm$meta$pred_spec$P_random > 0
}

has_random_scale <- function(lmmelsm) {
    lmmelsm$meta$pred_spec$Q_random > 0
}

get_number_re <- function(lmmelsm) {
    F <- get_F(lmmelsm)
    P_random <- get_P_random(lmmelsm)
    Q_random <- get_Q_random(lmmelsm)

    2 * F + P_random * F + Q_random * F
}

get_re_indices <- function(lmmelsm) {
    F <- lmmelsm$meta$indicator_spec$F
    P_random <- lmmelsm$meta$pred_spec$P_random
    Q_random <- lmmelsm$meta$pred_spec$Q_random

    mu_random <- 1:F
    logsd_random <- (F+1):(2*F)
    mu_beta_random <- if(has_random_location(lmmelsm)) {
                          (2*F + 1):(2*F + P_random*F)
                      } else {NA}
    logsd_beta_random <- if(has_random_scale(lmmelsm)) {
                             (2*F + P_random*F + 1):(2*F + P_random*F + Q_random*F)
                         } else {NA}

    nlist(mu_random,
          logsd_random,
          mu_beta_random,
          logsd_beta_random)
}

get_factor_names <- function(lmmelsm) {
    unlist(lmmelsm$meta$indicator_spec$fname)
}

get_indicator_names <- function(lmmelsm) {
    lmmelsm$meta$indicator_spec$mname
}

get_predictor_names <- function(lmmelsm, which = c("location", "scale", "between")) {
    mod <- match.arg(which)
    lmmelsm$meta$pred_spec$pname[[mod]]
}

get_group_name <- function(lmmelsm) {
    lmmelsm$meta$group_spec$name
}

get_group_numeric <- function(lmmelsm) {
    lmmelsm$meta$group_spec$numeric
}

get_group_labels <- function(lmmelsm) {
    lmmelsm$meta$group_spec$data
}

get_group_map <- function(lmmelsm) {
    lmmelsm$meta$group_spec$map
}

get_K <- function(lmmelsm) {
    lmmelsm$meta$group_spec$K
}

get_F <- function(lmmelsm) {
    lmmelsm$meta$indicator_spec$F
}

get_J <- function(lmmelsm) {
    lmmelsm$meta$indicator_spec$J
}

get_P <- function(lmmelsm) {
    lmmelsm$meta$pred_spec$P
}

get_Q <- function(lmmelsm) {
    lmmelsm$meta$pred_spec$Q
}

get_R <- function(lmmelsm) {
    lmmelsm$meta$pred_spec$R
}

get_P_random <- function(lmmelsm) {
    lmmelsm$meta$pred_spec$P_random
}

get_Q_random <- function(lmmelsm) {
    lmmelsm$meta$pred_spec$Q_random
}

get_P_random_ind <- function(lmmelsm) {
    lmmelsm$meta$pred_spec$P_random_ind
}

get_Q_random_ind <- function(lmmelsm) {
    lmmelsm$meta$pred_spec$Q_random_ind
}

get_N <- function(lmmelsm) {
    lmmelsm$meta$indicator_spec$N
}

has_multivariate <- function(lmmelsm) {
    lmmelsm$meta$indicator_spec$F > 1
}

get_S <- function(lmmelsm) {
    (lmmelsm$fit@stan_args[[1]]$iter - lmmelsm$fit@stan_args[[1]]$warmup) * get_number_chains(lmmelsm)
}

get_number_chains <- function(lmmelsm) {
    lmmelsm$stan_args$chains
}
