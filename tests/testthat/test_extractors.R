test_that("Object extractors work properly", {
    library(LMMELSM)
    data(sim_data)

    fit_obs <- lmmelsm(observed ~ A_1 + A_2, subject, sim_data, iter = 2, cores = 1, chains = 1)

    expect_equal(has_latent(fit_obs), FALSE)
    expect_equal(has_location(fit_obs), FALSE)
    expect_equal(has_scale(fit_obs), FALSE)
    expect_equal(has_between(fit_obs), FALSE)
    expect_equal(has_random_location(fit_obs), FALSE)
    expect_equal(has_random_scale(fit_obs), FALSE)
    expect_equal(get_number_re(fit_obs), 4)
    expect_equal(get_re_indices(fit_obs), list(mu_random = 1:2, logsd_random = 3:4, mu_beta_random = NA, logsd_beta_random = NA))
    expect_equal(has_multivariate(fit_obs), TRUE)

    fit_lat <- lmmelsm(list(A ~ A_1 + A_2 + A_3,
                            B ~ N_1 + N_2 + N_3,
                            location ~ x1,
                            scale ~ x2),
                       subject, sim_data, iter = 2, cores = 1, chains = 1)

    expect_equal(has_latent(fit_lat), TRUE)
    expect_equal(has_location(fit_lat), TRUE)
    expect_equal(has_scale(fit_lat), TRUE)
    expect_equal(has_between(fit_lat), FALSE)
    expect_equal(has_random_location(fit_lat), FALSE)
    expect_equal(has_random_scale(fit_lat), FALSE)
    expect_equal(get_number_re(fit_lat), 4)
    expect_equal(get_re_indices(fit_lat), list(mu_random = 1:2, logsd_random = 3:4, mu_beta_random = NA, logsd_beta_random = NA))
    expect_equal(has_multivariate(fit_lat), TRUE)

    fit_lat <- lmmelsm(list(A ~ A_1 + A_2 + A_3,
                            B ~ N_1 + N_2 + N_3,
                            location ~ x1 | x1,
                            scale ~ x1 + x2 | x1 + x2,
                            between ~ baseline),
                       subject, sim_data, iter = 2, cores = 1, chains = 1)

    expect_equal(has_latent(fit_lat), TRUE)
    expect_equal(has_location(fit_lat), TRUE)
    expect_equal(has_scale(fit_lat), TRUE)
    expect_equal(has_between(fit_lat), TRUE)
    expect_equal(has_random_location(fit_lat), TRUE)
    expect_equal(has_random_scale(fit_lat), TRUE)
    expect_equal(get_number_re(fit_lat), 4 + 2 + 4)
    expect_equal(get_re_indices(fit_lat), list(mu_random = 1:2, logsd_random = 3:4, mu_beta_random = 5:6, logsd_beta_random = 7:10))
    expect_equal(has_multivariate(fit_lat), TRUE)

    fit_obs <- lmmelsm(observed ~ A_1 + A_2, subject, sim_data, iter = 10, cores = 1, chains = 2)
    expect_equal(get_S(fit_obs), 10)
    expect_equal(get_N(fit_obs), nrow(sim_data))
    expect_equal(get_K(fit_obs), length(unique(sim_data$subject)))
})
