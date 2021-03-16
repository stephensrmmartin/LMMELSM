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


test_that("predict.lmmelsm is returning correct structures", {
    library(LMMELSM)
    data(sim_data)

    new_data <- sim_data[c(1,51,101),]

    iter <- 5
    cores <- 1
    chains <- 1

    fit_int <- lmmelsm(list(A ~ A_1 + A_2 + A_3 + A_4 + A_5 + A_6,
                            N ~ N_1 + N_2 + N_3 + N_4 + N_5 + N_6
                            ),
                       subject, sim_data, iter = iter, cores = cores, chains = chains)
    pred <- predict(fit_int, new_data)

    expect_equal(length(pred), 3)
    expect_equal(nrow(pred$eta), 3 * 2)
    expect_equal(nrow(pred$eta_logsd), 3 * 2)
    expect_equal(nrow(pred$y), 3 * 12)

    fit_obs <- lmmelsm(list(observed ~ A_1 + N_1), subject, sim_data, iter = iter, cores = cores, chains = chains)

    pred <- predict(fit_obs, new_data)

    expect_equal(nrow(pred$eta), nrow(pred$y))
    expect_equal(nrow(pred$eta), 3 * 2)
    expect_equal(nrow(pred$eta_logsd), 3 * 2)
    expect_equal(nrow(pred$y), 3 * 2)

    fit_obs_1 <- lmmelsm(observed ~ A_1, subject, sim_data, iter = iter, cores = cores, chains = chains)

    pred <- predict(fit_obs_1, new_data)

    expect_equal(nrow(pred$eta), nrow(pred$y))
    expect_equal(nrow(pred$eta), 3 * 1)
    expect_equal(nrow(pred$eta_logsd), 3 * 1)
    expect_equal(nrow(pred$y), 3 * 1)

    # Testing full gambit
    fit_lat <- lmmelsm(list(A ~ A_1 + A_2 + A_3 + A_4 + A_5 + A_6,
                            N ~ N_1 + N_2 + N_3 + N_4 + N_5 + N_6,
                            location ~ x1 + x2 + baseline | x1 + x2,
                            scale ~ x1 + x2 + baseline | x1 + x2,
                            between ~ baseline),
                       subject, sim_data, iter = iter, cores = cores, chains = chains)

    pred <- predict(fit_lat, new_data)

    expect_equal(length(pred), 3)
    expect_equal(nrow(pred$eta), 3 * 2)
    expect_equal(nrow(pred$eta_logsd), 3 * 2)
    expect_equal(nrow(pred$y), 3 * 12)

    pred <- predict(fit_lat, new_data, summarize = FALSE)
    expect_equal(length(pred), 3)
    expect_equal(length(pred[[1]]), 3)
    expect_equal(ncol(pred[[1]][["eta"]]), 2)
    expect_equal(ncol(pred[[1]][["eta_logsd"]]), 2)
    expect_equal(ncol(pred[[1]][["y"]]), 12)
})

test_that("fitted returns correct structures", {
    library(LMMELSM)
    data(sim_data)

    iter <- 5
    cores <- 1
    chains <- 1

    fit_lat <- lmmelsm(list(A ~ A_1 + A_2 + A_3 + A_4 + A_5 + A_6,
                            N ~ N_1 + N_2 + N_3 + N_4 + N_5 + N_6,
                            location ~ x1 + x2 + baseline | x1 + x2,
                            scale ~ x1 + x2 + baseline | x1 + x2,
                            between ~ baseline),
                       subject, sim_data, iter = iter, cores = cores, chains = chains)

    fitted_lat <- fitted(fit_lat)
    expect_equal(length(fitted_lat), 2)
    expect_equal(nrow(fitted_lat[[1]]), nrow(sim_data) * 2)
    expect_equal(nrow(fitted_lat[[2]]), nrow(sim_data) * 2)
})
