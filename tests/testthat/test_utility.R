test_that(".list_zip", {
    test_list <- list(A = list(1,2,3), B = list(4,5,6), C = list(7,8,9))
    list_zip_output <- list(t(c(1,4,7)), t(c(2,5,8)), t(c(3,6,9)))
    list_zip_output <- lapply(list_zip_output, function(l) {colnames(l) <- c("A","B","C"); l})
    expect_equal(do.call(.list_zip, test_list), list_zip_output)
})
