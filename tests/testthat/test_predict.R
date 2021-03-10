test_that(".add_group_codings returns numerics and NAs.", {
    facLabels <- c("a","b","c","d","e")
    gs <- list(name = "subject", data = as.factor(facLabels),
               numeric = as.numeric(as.factor(facLabels)),
               K = 5, map = data.frame(numeric = as.numeric(as.factor(facLabels)),
                                       label = as.factor(facLabels)))

    newdata <- data.frame(subject = c("c", "f", NA),
                          A_1 = rnorm(3))

    newdata <- .add_group_codings(gs, newdata)
    expect_equal(newdata[1, "subject"], 3)
    expect_equal(is.na(newdata[2, "subject"]), TRUE)
    expect_equal(is.na(newdata[3, "subject"]), TRUE)

    newdata <- data.frame(A_1 = rnorm(3))
    newdata <- .add_group_codings(gs, newdata)
    for(n in 1:nrow(newdata)) {
        expect_equal(is.na(newdata[n, "subject"]), TRUE)
    }
})
