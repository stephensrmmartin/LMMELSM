test_that("lmmelsm returns correct structures.", {
    set.seed(13)
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

    d <- LMMELSM:::simulate_lmmelsm(
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

    fit <- lmmelsm(list(fac1 ~ obs_1 + obs_2 + obs_3 + obs_4,
                        fac2 ~ obs_5 + obs_6 + obs_7 + obs_8,
                        location ~ loc_1 + loc_2 | loc_1,
                        scale ~ sca_1 + sca_2 | sca_2),
                group = subject, data = d$df, iter = 10, chains = 1, cores = 1)

    expect_s3_class(fit, "lmmelsm")
    expect_named(fit, c("fit", "meta", "data","stan_data", "stan_args"))

    expect_equal(fit$meta$group_spec$K, K)
    expect_equal(fit$meta$indicator_spec$J, J)
    expect_equal(fit$meta$indicator_spec$N, n*K)

    fit_sum <- summary(fit)
    expect_s3_class(fit_sum, "summary.lmmelsm")
    expect_named(fit_sum, c("meta", "summary"))
    expect_named(fit_sum$summary,
		c("lambda",
		"nu",
		"sigma",
                "mu_coef",
                "logsd_coef",
                "zeta",
		"random_mu_intercept",
		"random_logsd_intercept",
		"random_mu_coef",
		"random_logsd_coef",
		"random_sigma",
		"random_correlation",
		"factor_correlation"),
		ignore.order = TRUE)

    fit_ranef <- ranef(fit)
    expect_type(fit_ranef, "list")
    expect_named(fit_ranef, c("random_mu_intercept",
                              "random_logsd_intercept",
                              "random_mu_coef",
                              "random_logsd_coef"), ignore.order = TRUE)
    expect_equal(nrow(fit_ranef$random_mu_intercept), K*F)
    expect_equal(nrow(fit_ranef$random_logsd_intercept), K*F)

    fit_coef <- coef(fit)
    expect_named(fit_coef, c("mu_intercept",
                             "logsd_intercept",
                             "mu_coef",
                             "logsd_coef"))
    expect_equal(nrow(fit_coef$mu_intercept), K*F)
    expect_equal(nrow(fit_coef$logsd_intercept), K*F)
})
