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
##' @title 
##' @param x lmmelsm object.
##' @param ... Not used.
##' @return NULL.
##' @author Stephen R. Martin
##' @export
print.lmmelsm <- function(x, ...) {
    .sep()

    cat("Estimated on:")
    cat(x$fit@date)
    .newline()
    cat("Elapsed time (seconds): ")
    cat(.get_elapsed_time(x))

    .sep()

    cat("Model specification:")
    .newline(2)

    cat("Measurement model")
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

summary.lmmelsm <- function(object, ...) {
    
}

print.summary.lmmelsm <- function(x, ...) {
    
}

.print.formula <- function(f) {
    string <- Reduce(paste0, deparse(f))
    return(string)
}

.get_elapsed_time <- function(sOut) {
    times <- rstan::get_elapsed_time(sOut$fit)
    out <- max(rowSums(times))

    return(out)
}
