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
    if(!is.matrix(s) & inherits(s, 'stanfit')) {
        s <- as.matrix(s, pars = par)
    }
    cns <- colnames(s)
    S <- nrow(s)

    ## rexInner <- r"(.*\[(\d+(?:,\d+)*)\])" # R < 4.0 breaks
    rexInner <-  ".*\\[(\\d+(?:,\\d+)*)\\]"
    inner <- gsub(rexInner, "\\1", cns)

    inds <- apply(do.call(rbind, strsplit(inner, ",")), 2, as.numeric)
    if(is.null(dim(inds))) { # Edge case; only one set of indices
        dim(inds) <- c(1,length(inds))
    }
    inds_max <- apply(inds, 2, max)
    inds_all <- c(inds_max, S)

    arr <- array(t(s), dim = inds_all)
    return(arr)
}

# Turns [,,..., S] array into S-length list of [,,...]
.array_to_list <- function(arr) {
    dims <- dim(arr)
    all_dims <- lapply(dims[-length(dims)], function(d){1:d})
    S <- dims[length(dims)]

    lst <- lapply(seq_len(S), function(s) {
        args <- c(list(arr), all_dims, s, drop = FALSE)
        do.call(`[`, args)
    })

    lst <- lapply(lst, function(s){
        array(s, dims[-length(dims)])
    })

    lst
}

.extract_transform_to_list <- function(s, par = NULL) {
    # Get as array[,,...,S]
    arr <- .extract_transform(s, par)
    lst <- .array_to_list(arr)

    lst
}

.ones <- function(x) {
    out <- matrix(1, x, 1)
    return(out)
}
##' @title Zip two lists together with function.
##' @param f Function
##' @param ... Lists to zip.
##' @return List.
##' @author Stephen Martin
##' @keywords internal
.list_zip <- function(..., f = cbind) {
    dots <- list(...)
    lengths <- lapply(dots, length)
    ## stopifnot(all.equal(lengths))

    lapply(seq_len(lengths[[1]]), function(i) {
        items <- lapply(dots, function(l) {
            l[[i]]
        })
        do.call(f, items)
    })
}
