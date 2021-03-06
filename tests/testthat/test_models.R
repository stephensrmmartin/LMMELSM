test_that("Specs return correct model dimensions", {
    library(LMMELSM)
    data(sim_data)

    df <- sim_data
    df <- df[df$subject %in% 1:20,]

    unidim <- list(A ~ A_1 + A_2 + A_3 + A_4 + A_5 + A_6)
    multidim <- list(A ~ A_1 + A_2 + A_3 + A_4 + A_5 + A_6,
                     B ~ N_1 + N_2 + N_3 + N_4 + N_5 + N_6)

    unidim_X1 <- c(unidim, location ~ baseline)
    multidim_X1 <- c(multidim, location ~ baseline)

    unidim_obs <- list(observed ~ A_1)
    multidim_obs <- list(observed ~ A_1 + A_2 + A_3 + A_4 + A_5 + A_6)

    ## unidim__X1 <- c(unidim, scale ~ baseline)
    ## multidim__X1 <- c(multidim, scale ~ baseline)

    ## unidim_x1 <- c(unidim, location ~ x1)
    ## multidim_x1 <- c(multidim, location ~ x1)

    ## unidim__x1 <- c(unidim, scale ~ x1)
    ## multidim__x1 <- c(multidim, scale ~ x1)

    ## unidim_x1re <- c(unidim, location ~ x1|x1)
    ## multidim_x1re <- c(multidim, location ~ x1|x1)

    ## unidim__x1re <- c(unidim, scale ~ x1|x1)
    ## multidim__x1re <- c(multidim, scale ~ x1|x1)

    unidim_full <- c(unidim, location ~ x1 + baseline|x1, scale ~ x2 + baseline|x2, between ~ baseline)
    multidim_full <- c(multidim, location ~ x1 + baseline|x1, scale ~ x2 + baseline|x2, between ~ baseline)

    form_names <- ls(pattern = "(uni|multi)dim.*")
    forms <- mget(form_names)

    iter <- 5
    cores <- 1
    chains <- 1

    fits <- lapply(forms, function(x) {
        lmmelsm(x, subject, df, iter = iter, cores = cores, chains = chains)
    })
    for(f in fits) {
        expect_s3_class(f, "lmmelsm")
    }

    ###################
    # Tests           #
    ###################
    for(f in fits) {
        f_sum <- summary(f)
        f_ranef <- ranef(f)
        f_coef <- coef(f)

        # Summary class
        expect_s3_class(f_sum, "summary.lmmelsm")

        # L2_pred_only should be TRUE if x_loc and x_sca only have (less than) K unique values; except in bizarre edge cases.
        expect_equal(f_sum$meta$pred_spec$L2_pred_only, length(unique(f$meta$pred_spec$x_loc)) <= f$meta$group_spec$K & length(unique(f$meta$pred_spec$x_sca)) <= f$meta$group_spec$K)

        # ranef intercepts should have F * K rows.
        expect_equal(nrow(f_ranef$random_mu_intercept), f$meta$group_spec$K * f$meta$indicator_spec$F)
        expect_equal(nrow(f_ranef$random_logsd_intercept), f$meta$group_spec$K * f$meta$indicator_spec$F)

        # coef intercepts should also have F * K rows.
        expect_equal(nrow(f_coef$mu_intercept), f$meta$group_spec$K * f$meta$indicator_spec$F)
        expect_equal(nrow(f_coef$logsd_intercept), f$meta$group_spec$K * f$meta$indicator_spec$F)

        # ranef and coefs should be 4-length lists
        expect_length(f_coef, 4)
        expect_length(f_ranef, 4)

        # coefs should all equal ranefs for the intercepts
        for(i in seq_len(f$meta$indicator_spec$F * f$meta$group_spec$K)) {
            if(f$meta$latent) {
                expect_equal(f_ranef$random_mu_intercept[i,"Mean"], f_coef$mu_intercept[i,"Mean"])
                expect_equal(f_ranef$random_logsd_intercept[i,"Mean"], f_coef$logsd_intercept[i,"Mean"])
            }
            if(!f$meta$latent) { # And equal once subtracting off fixed effects
                expect_equal(f_ranef$random_mu_intercept[i, "Mean"], f_coef$mu_intercept[i, "Mean"] - rep(f_sum$summary$nu[,"Mean"], each = f$meta$group_spec$K)[i])
            }
        }
    }

})
