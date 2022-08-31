##' @title Multiply a row by a list of matrices
##' @param r Row vector.
##' @param mats 
##' @return List of row by matrix products.
##' @author Stephen R Martin
##' @keywords internal
row_multiply_list_mats <- function(r, mats) {
    r <- matrix(r, nrow = 1)
    lapply(mats, function(m) {
        r %*% m
    })
}

list_extract <- function(l, what) {
    lapply(l, function(l) {
        l[[what]]
    })
}
