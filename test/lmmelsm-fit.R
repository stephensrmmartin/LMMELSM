library(LMMELSM)

data(sim_data)

sub_n <- 20

get_subsample_indices <- function(group, sub_n) {
    indices <- seq_len(length(group))

    indices_by_group <- tapply(indices, group, function(x) {
        sample(x, sub_n)
    })

    do.call(c, indices_by_group)
}

set.seed(13)

indices <- get_subsample_indices(sim_data$subject, sub_n)
ds <- sim_data[indices, ]

# Univariate #

fit <- lmmelsm(list(A ~ A_1 + A_2 + A_3 + A_4 + A_5 + A_6), subject, ds, iter = 1000)

# Univariate_observed #
ds$A_mean <- rowMeans(ds[,paste0("A_",1:6)])
fit_uni_obs <- lmmelsm(observed ~ A_mean, subject, ds, iter = 1000)

# Multivariate_observed #

fit_multi_obs <- lmmelsm(observed ~ A_1 + A_2 + A_3 + A_4 + A_5 + A_6, subject, ds, iter = 1000)

# Multivariate #

fit_multi <- lmmelsm(list(A ~ A_1 + A_2 + A_3 + A_4 + A_5 + A_6,
                          N ~ N_1 + N_2 + N_3 + N_4 + N_5 + N_6),
                     subject, ds, iter = 1000)



# True model #
fit_true <- lmmelsm(list(A ~ A_1 + A_2 + A_3 + A_4 + A_5 + A_6,
                         N ~ N_1 + N_2 + N_3 + N_4 + N_5 + N_6,
                         location ~ x1 + baseline | x1,
                         scale ~ x2 + baseline | x2,
                         between ~ baseline),
                    subject, ds, iter = 1000)


# Predict #
pred_fit_true <- predict(fit_true, include_error = FALSE)
fit_true_fit_true <- fitted(fit_true)

head(d$params$eta)
head(d$params$eta[indices,])
