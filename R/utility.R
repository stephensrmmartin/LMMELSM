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
