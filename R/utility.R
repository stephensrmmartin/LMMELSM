##' @title Operator for testing NULL and returning expr if NULL
##' @param object Object to test for NULL
##' @param expr Expression to evaluate and return if object is NULL
##' @return object if object is non-NULL, expression output if object is NULL.
##' @author Stephen R. Martin
##' @keywords internal
`%IfNull%` <- function(object, expr) {
    if(is.null(object)) {
        return(eval(expr))
    } else {
        return(object)
    }
}

##' @title Creates named list.
##' @param ... Objects for list.
##' @return Named List.
##' @author Stephen R. Martin
##' @keywords internal
nlist <- function(...) {
    mc <- match.call()
    out <- list(...)

    not_named <- is.null(names(out))
    is_named <- if(not_named) {
                    FALSE
                } else {
                    nzchar(names(out))
                }

    args <- as.character(mc)[-1] # Not the fn name.

    if(not_named) {
        names(out) <- args
    } else {
        names(out)[!is_named] <- args[!is_named]
    }

    return(out)
}

.array_extract <- function(a, ind) {
    dim_a <- dim(a)
    lastDim <- length(dim_a)

    dims <- lapply(dim_a[1:(lastDim - 1)], function(x) {
        1:x
    })

    args <- c(list(a), dims, ind, drop = FALSE)
    a_sub <- do.call(`[`, args)

    out <- array(a_sub, dim = dim_a[1:(lastDim - 1)])

    return(out)
}

.prob_to_interval <- function(p) {
    low <- (1 - p) / 2
    high <- 1 - low
    return(c(low, high))
}

.extract_transform <- function(s, par = NULL) {
    if(!is.matrix(s) & class(s) == "stanfit") {
        s <- as.matrix(s, pars = par)
    }
    cns <- colnames(s)
    S <- nrow(s)

    ## rexInner <- r"(.*\[(\d+(?:,\d+)*)\])" # R < 4.0 breaks
    rexInner <-  ".*\\[(\\d+(?:,\\d+)*)\\]"
    inner <- gsub(rexInner, "\\1", cns)

    inds <- apply(do.call(rbind, strsplit(inner, ",")), 2, as.numeric)
    inds_max <- apply(inds, 2, max)
    inds_all <- c(inds_max, S)

    arr <- array(t(s), dim = inds_all)
    return(arr)
}

.ones <- function(x) {
    out <- matrix(1, x, 1)
    return(out)
}
