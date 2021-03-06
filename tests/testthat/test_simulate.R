test_that("Simulated unidimensional data are correct.", {
    d <- simulate_lmmelsm(
        n = 50,
        K = 10,
        lambda = c(.7,.7,.7,.8,.9),
        resid = rep(1, 5),
        nu = rep(0, 5),
        mu_beta = matrix(c(.4, -.6), ncol = 1),
        logsd_beta = matrix(c(.4, -.6), ncol = 1),
        P_random_ind = 1:2,
        Q_random_ind = 1:2,
        mu_logsd_betas_cor = diag(1, 2 + 2 + 2),
        mu_logsd_betas_sigma = rep(.3, 2 + 2 + 2),
        epsilon_cor = matrix(1, 1, 1)
    )

    expect_named(d, c("params", "data", "df"))
    expect_s3_class(d$df, "data.frame")
    df <- d$df
    expect_equal(nrow(df), d$params$N)
    expect_named(df, c(paste0("obs_",1:5), "subject", paste0("loc_",1:2), paste0("sca_", 1:2)))
    expect_equal(d$params$lambda, matrix(c(.7,.7,.7,.8,.9), nrow = 1))
    print(names(d$params))
    expect_equal(nrow(d$params$mu_logsd_betas_re), 10)
    expect_equal(ncol(d$params$mu_logsd_betas_re), 2 * 1 + 2 * 1 + 2 * 1)
})


test_that("Simulated multidimensional data are correct.", {
    n <- 50
    K <- 10
    F <- 2
    J <- 8
    lambda <- matrix(c(.8,.8,.8,.8,0,0,0,0,
                       0,0,0,0,.8,.8,.8,.8),byrow=TRUE,nrow = F)
    resid <- rep(1, J)
    nu <- rep(0, J)
    mu_beta <- matrix(c(.4,.5,.6,.7), ncol = F)
    logsd_beta <- matrix(c(.4,.5,.6,.7), ncol = F)
    P_random_ind <- 1
    Q_random_ind <- 2
    mu_logsd_betas_cor <- diag(1, 2 * F + F + F)
    mu_logsd_betas_sigma <- rep(.3, 2 * F + F + F)
    epsilon_cor <- diag(1,2,2)
    d <- simulate_lmmelsm(
        n = n,
        K = K,
        lambda = lambda,
        resid = resid,
        nu = nu,
        mu_beta = mu_beta,
        logsd_beta = logsd_beta,
        P_random_ind = P_random_ind,
        Q_random_ind = Q_random_ind,
        mu_logsd_betas_cor = mu_logsd_betas_cor,
        mu_logsd_betas_sigma = mu_logsd_betas_sigma,
        epsilon_cor = epsilon_cor
    )

    expect_equal(d$params$N, n*K)
    expect_equal(d$params$P, 2)
    expect_equal(d$params$Q, 2)
    expect_equal(d$params$P_random, 1)
    expect_equal(d$params$Q_random, 1)
    expect_equal(d$params$P_random_ind, 1)
    expect_equal(d$params$Q_random_ind, 2)
    expect_equal(ncol(d$params$eta), 2)
    expect_equal(ncol(d$params$mu_logsd_betas_re), F + F + 1*F + 1*F)

    df <- d$df
    expect_equal(nrow(df), n*K)
    expect_equal(ncol(df), J + 1 + 2 + 2)
})
