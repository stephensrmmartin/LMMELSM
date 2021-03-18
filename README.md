
# LMMELSM

<!-- badges: start -->
<!-- badges: end -->

LMMELSM stands for latent, multidimensional mixed effects location scale models.

The mixed effects location scale model is a multilevel model in which both the location ($\mu$) and scale ($\sigma$) are modeled and have random effects.
That is, in *addition* to random intercepts and coefficients, the MELSM allows variance to differ between groups, and across covariates.
Within the MELSM framework, variances are meaningful, and can be directly modeled.
MELSMs are useful for predicting volatility, variance, uncertainty, variation between group means, measurement error variance, and much more.
As an added side effect, prediction intervals are tailored to each individual or group (differing predictive uncertainty), and fixed effects are more robust because those with less variance (more reliability) are implicitly weighted more.

LMMELSM extends the MELSM to include *latent* variables.
Instead of predicting the expected values and variances of observed or computed variables, the LMMELSM predicts expected values and variances of multiple correlated latent factors.
Assuming an adequate measurement model, the LMMELSM therefore models the variance of a reduced-error score.

LMMELSM is therefore useful for modeling conditional and group-specific latent means and variances.
One such example is modeling the intraindividual variability of latent scores of a persons over time in an experience sampling methodology.
Another example is modeling the heterogeneity of individuals' latent means from person-level covariates.

LMMELSM currently supports:
- Observed outcomes
- Latent outcomes
- Single or multiple outcomes
- Location modeling, with random intercepts and coefficients
- Within-group scale modeling, with random intercepts and coefficients
- Between-group scale modeling (on the between-group mean variances)

## Installation

You can install the released version of LMMELSM from [github](https://github.com/stephensrmmartin/LMMELSM) with:

``` r
remotes::install_github("stephensrmmartin/LMMELSM")
```

## Example

In this example, increasingly complex MELSMs are modeled, and demonstrate the simple formula syntax.
``` r
library(LMMELSM)

data(sim_data)

# Intercept-only model
fit <- lmmelsm(Agreeableness ~ A_1 + A_2 + A_3 + A_4 + A_5 + A_6, subject, sim_data)

# Random Time-varying predictor on location
fit <- lmmelsm(list(Agreeableness ~ A_1 + A_2 + A_3 + A_4 + A_5 + A_6,
                    location ~ x1 | x1),
               subject, sim_data)

# Random Time-varying predictors on both location and scale
fit <- lmmelsm(list(Agreeableness ~ A_1 + A_2 + A_3 + A_4 + A_5 + A_6,
                    location ~ x1 | x1,
                    scale ~ x2 | x2),
               subject, sim_data)

# Time-varying predictors, person-level predictors on location, scale, and between-group variance
# Multidimensional!
fit <- lmmelsm(list(Agreeableness ~ A_1 + A_2 + A_3 + A_4 + A_5 + A_6,
                    Neuroticism ~ N_1 + N_2 + N_3 + N_4 + N_5 + N_6,
                    location ~ x1 + baseline | x1,
                    scale ~ x2 + baseline | x2,
                    between ~ baseline),
               subject, sim_data)

# Time-varying predictors, person-level predictors on location, scale, and between-group variance
fit <- lmmelsm(list(Agreeableness ~ A_1 + A_2 + A_3 + A_4 + A_5 + A_6,
                    Neuroticism ~ N_1 + N_2 + N_3 + N_4 + N_5 + N_6,
                    location ~ x1 + baseline | x1,
                    scale ~ x1 + x2 + baseline | x1 + x2,
                    between ~ baseline),
               subject, sim_data)


# Non-latent Multivariate MELSM
fit <- lmmelsm(list(observed ~ A_1 + N_1,
                    location ~ x1 + baseline|x1,
                    scale ~ x2 + baseline|x2,
                    between ~ baseline),
               subject, sim_data)

# Summarize
summary(fit)

# Get subject-specific random effects
ranef(fit)

# Get subject-specific coefficients
coef(fit)

```

