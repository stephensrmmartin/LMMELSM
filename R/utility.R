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
