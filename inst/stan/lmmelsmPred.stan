functions {

  // Takes indicator spec and lambda_est, unrolls into lambda matrix.
  matrix lambda_mat(int J, int F, array[] int J_f, array[,] int F_ind, vector lambda_est) {
    matrix[F, J] out = rep_matrix(0.0, F, J);
    int count = 1;
    for(f in 1:F) {
      for(jj in 1:J_f[f]) {
	out[f, F_ind[f, jj]] = lambda_est[count];
	count += 1;
	}
    }
    return(out);
  }

  /*
    Non-centered RE params to REs
    @param matrix[K, num] z; K = groups, num = number of REs
    @param cholesky_factor_corr[num] L;
    @param vector[num] sigma; RE SDs
    @return matrix[K, num] REs;
   */
  matrix z_to_re(matrix z, matrix L, vector sigma) {
    return(z * diag_pre_multiply(sigma, L)');
  }

  /*
    Non-centered RE params to REs, with predicted variances.
    @param matrix[K, num] z; K = groups, num = number of REs.
    @param cholesk_factor_corr[num] L;
    @param vector[num] sigma; RE SDs (intercept, on exp scale)
    @param matrix[K, R] x_bet_l2; Level-2 design matrix for log-scale contributions to bet. group SDs.
    @param matrix[R, num] zeta; coefficient matrix.
    @param matrix[K, num] REs;
   */
  matrix z_to_re_bet(matrix z, matrix L, vector sigma, matrix x_bet_l2, matrix zeta) {
    int num = cols(z);
    int K = rows(z);
    matrix[K, num] sigmas = (rep_vector(1.0, K) * sigma') .* exp(x_bet_l2 * zeta);
    matrix[K, num] out;

    for(k in 1:K) {
      out[k] = z[k] * diag_pre_multiply(sigmas[k], L)';
    }

    return(out);
  }

  /*
    Non-centered RE params to REs, with predicted variances for the intercepts.
    @param matrix[K, num] z; K = groups, num = number of REs.
    @param cholesk_factor_corr[num] L;
    @param vector[num] sigma; RE SDs (intercept, on exp scale)
    @param matrix[K, R] x_bet_l2; Level-2 design matrix for log-scale contributions to bet. group SDs.
    @param matrix[R, re_intercepts] zeta; coefficient matrix.
    @param matrix[K, num] REs;
   */
  matrix z_to_re_bet_intercepts(matrix z, matrix L, vector sigma, matrix x_bet_l2, matrix zeta) {
    int re_total = cols(z);
    int K = rows(z);
    int re_intercepts = cols(zeta);
    matrix[K, re_total] sigmas = (rep_vector(1.0, K) * sigma');
    matrix[K, re_total] out;

    sigmas[, 1:re_intercepts] .*= exp(x_bet_l2 * zeta);

    for(k in 1:K) {
      out[k] = z[k] * diag_pre_multiply(sigmas[k], L)';
    }

    return(out);
  }

  /*
    Convert repeated measures to subject-level dataset
    @return matrix[K, cols(l1)]; in order of 1:K, if using l1_to_l2_indices().
   */
  matrix l1_to_l2(matrix l1, array[] int indices) {
    int K = size(indices);
    int n_col = cols(l1);
    matrix[K, n_col] l2 = l1[indices];
    return(l2);
  }

  /*
    Find indices for which each group, k in 1:K, first appears in array[] int group.
    @param int K: Total number of groups.
    @param array[] int group: Array of group integers from L1 dataset.
    @return array[K] int; Array of L1 indices in which each k first appears.
   */
  array[] int l1_to_l2_indices(int K, array[] int group) {
    int N = size(group);
    array[K] int where_l1_first_k = rep_array(0, K);

    for(n in 1:N) {
      if(where_l1_first_k[group[n]] == 0) {
	where_l1_first_k[group[n]] = n;
      }
    }

    return(where_l1_first_k);
  }

  /*
    Given a KxM matrix, convert to a K-length array of RxC matrices.
    @param int R: Rows in output matrices.
    @param int C: Columns in output matrices.
    @param matrix mat: KxM Matrix to restructure.
    @return array[] matrix: Array of RxC matrices.
   */
  array[] matrix mat_to_mat_array(int R, int C, matrix mat) {
    int K = rows(mat);
    array[K] matrix[R, C] out;

    for(k in 1:K) {
      out[k] = to_matrix(mat[k], R, C);
    }

    return(out);
  }

  /*
    Creates integer sequence. Useful for pre-defining indices.
    Conceptually, same as 'from:to' in R; but Stan does not have this.
    @param int from
    @param int to
    @return array[] int: Array of integers in sequence defined by from, to; inclusive.
   */
  array[] int seq_from_to(int from, int to) {
    int length = to - from + 1;
    array[length] int out;
    for(i in 1:length) {
      out[i] = from + i - 1;
    }

    return(out);
  }

}

data {
  int N; // Number of observations
  int J; // Number of indicators
  int F; // Number of latent factors
  int K; // Number of groups
  int P; // Number of predictors (location); [Does not include intercept]
  int Q; // Number of predictors (logsd); [Does not include intercept]
  int R; // Number of predictors (between-logsd); Fixed effects only (inherently a level-2 question)
  int P_random; // Number of random location coefficients
  int Q_random; // Number of random scale coefficients
  array[P_random] int P_random_ind; // Indices in x_loc corresponding to random location predictors
  array[Q_random] int Q_random_ind; // Indices in x_loc corresponding to random scale predictors

  array[N] int group; // Grouping indicator
  matrix[N, P] x_loc; // Location predictors
  matrix[N, Q] x_sca; // Scale predictors
  matrix[N, R] x_bet; // Between predictors

  // Indicators
  array[F] int J_f; // Number of indicators for each factor
  array[F, J] int F_ind; // Indicator indices.
  matrix[N, J] y; // Indicator data.

  // Options
  int<lower=0, upper=1> prior_only; // Whether to sample from prior only
  int<lower=0, upper=1> L2_pred_only; // Whether only level-2 predictors are provided; allows quicker computations.
  
}

transformed data {
  int intercept_only = P == 0 && Q == 0; // Whether intercept-only; allows quicker computations
  int lambda_total = sum(J_f);
  array[K] int l1_indices = l1_to_l2_indices(K, group);
  // Level 2 datasets, for ReVar and efficiency.
  matrix[K, P] x_loc_l2;
  matrix[K, Q] x_sca_l2;
  matrix[K, R] x_bet_l2;
  // Number of REs
  int re_intercepts = F*2;
  int re_mu_betas = P_random*F;
  int re_logsd_betas = Q_random*F;
  int re_total = re_intercepts + re_mu_betas + re_logsd_betas;
  // Pre-compute RE matrix indices
  array[F] int re_ind_mu = seq_from_to(1, F);
  array[F] int re_ind_logsd = seq_from_to(F + 1, F*2);
  array[re_mu_betas] int re_ind_mu_betas = seq_from_to((F*2) + 1, (F*2) + P_random*F);
  array[re_logsd_betas] int re_ind_logsd_betas = seq_from_to(F*2 + P_random*F + 1, F*2 + P_random*F + Q_random*F);
  // Pre-extract random design matrices.
  matrix[N, P_random] x_loc_re = x_loc[, P_random_ind];
  matrix[N, Q_random] x_sca_re = x_sca[, Q_random_ind];

  // Create L2 datasets for efficiency and ReVar model
  if(L2_pred_only) {
    x_loc_l2 = l1_to_l2(x_loc, l1_indices);
    x_sca_l2 = l1_to_l2(x_sca, l1_indices);
  }
  if(R > 0) {
    x_bet_l2 = l1_to_l2(x_bet, l1_indices);
  }

}

parameters {
  // Measurement model
  row_vector[J] nu;
  vector<lower=0>[lambda_total] lambda_est; // Assumes all loadings positive (including cross-loadings!)
  row_vector<lower=0>[J] sigma; // No error var scale model.

  // Structural model

  //// Location model: mu[ik] = B[i,0] + x[ik]B; B[i,0] = 0 + u[i,0]^(mu)
  matrix[P, F] mu_beta; // Fixed multivariate coefficients (Not including intercept)

  //// Scale model: eta[ik] = mu[ik] + epsilon[ik]; epsilon[ik] ~ mvn(0, f(L, logsd[ik])); logsd[ik] = G[i,0] + z[ik]G;
  matrix[N, F] epsilon_z; // Stochastic latent error
  matrix[Q, F] logsd_beta; // Fixed multivariate coefficients for scale model (Not including intercept)
  cholesky_factor_corr[F] epsilon_L; // Assumed equal across K.

  //// REs

  // For REs: We have P_random coefficients *per* factor; vector eta[ik] = X * matrix(B) + Z * matrix(u_i)
  matrix[K, re_total] mu_logsd_betas_random_z; // Random intercepts for mu, logsd
  cholesky_factor_corr[re_total] mu_logsd_betas_random_L;
  vector<lower=0>[re_total] mu_logsd_betas_random_sigma;

  // Between-group variance model
  // matrix[R, re_total] zeta;
  matrix[R, re_intercepts] zeta;
}

transformed parameters {
  matrix[F, J] lambda = lambda_mat(J, F, J_f, F_ind, lambda_est);
  matrix[K, re_total] mu_logsd_betas_random = R < 1 ?
    z_to_re(mu_logsd_betas_random_z, mu_logsd_betas_random_L, mu_logsd_betas_random_sigma) :
    z_to_re_bet_intercepts(mu_logsd_betas_random_z, mu_logsd_betas_random_L, mu_logsd_betas_random_sigma, x_bet_l2, zeta);
  matrix[K, F] mu_random = mu_logsd_betas_random[, re_ind_mu];
  matrix[K, F] logsd_random = mu_logsd_betas_random[, re_ind_logsd];
  array[K] matrix[P_random, F] mu_beta_random = mat_to_mat_array(P_random, F, mu_logsd_betas_random[, re_ind_mu_betas]);
  array[K] matrix[Q_random, F] logsd_beta_random = mat_to_mat_array(Q_random, F, mu_logsd_betas_random[,re_ind_logsd_betas]);
  matrix[N, F] eta;
  matrix[N, F] eta_logsd;
  // Location Predictions
  eta = mu_random[group]; // eta = 0 + u[0i]^(mu)
  if(P >= 1) { // + XB
    if(L2_pred_only) { // Multiply once, then broadcast; No random effects possible.
      eta += (x_loc_l2 * mu_beta)[group];
    } else {
      eta += x_loc * mu_beta;
    }
  }
  // Random effects
  if(P_random >= 1) {
    for(n in 1:N) {
      eta[n] += x_loc_re[n] * mu_beta_random[group[n]];
    }
  }
  // Scale Predictions
  eta_logsd = logsd_random[group]; // logsd = 0 + u[0i]^(logsd)
  if(Q >= 1) { // + ZG
    if(L2_pred_only) { // Multiply once, then broadcast; No random effects possible
      eta_logsd += (x_sca_l2 * logsd_beta)[group];
    } else {
      eta_logsd += x_sca * logsd_beta;
    }
  }
  // Random effects
  if(Q_random >= 1) {
    for(n in 1:N) {
      eta_logsd[n] += x_sca_re[n] * logsd_beta_random[group[n]];
    }
  }

  // Stochastic realizations
  if(L2_pred_only || intercept_only) { // Compute t(L_cov) for each k.
    {
      array[K] matrix[F, F] epsilon_cov_U;
      for(k in 1:K) {
	epsilon_cov_U[k] = (diag_pre_multiply(exp(eta_logsd[l1_indices[k]]), epsilon_L))';
      }
      for(n in 1:N) {
	eta[n] += epsilon_z[n] * epsilon_cov_U[group[n]];
      }

    }
  } else { // Compute t(L_cov) for each observation (much slower, but needed if sd varies by l1 predictor.)
    matrix[N, F] eta_sd = exp(eta_logsd);
    for(n in 1:N) {
      eta[n] += epsilon_z[n] * diag_pre_multiply(eta_sd[n], epsilon_L)';
    }
  }
  
}

model {

  // Priors
  nu ~ std_normal();
  lambda_est ~ std_normal();
  sigma ~ std_normal();

  to_vector(zeta) ~ std_normal();
  to_vector(mu_beta) ~ std_normal();
  to_vector(logsd_beta) ~ std_normal();
  to_vector(epsilon_z) ~ std_normal();
  epsilon_L ~ lkj_corr_cholesky(1);
  to_vector(mu_logsd_betas_random_z) ~ std_normal();
  mu_logsd_betas_random_L ~ lkj_corr_cholesky(1);
  mu_logsd_betas_random_sigma ~ std_normal();

  if(!prior_only){
    for(j in 1:J) {
      y[,j] ~ normal_id_glm(eta, nu[j], lambda[,j], sigma[j]);
    }
  }
}

generated quantities {
  corr_matrix[F] Omega_eta = multiply_lower_tri_self_transpose(epsilon_L);
  corr_matrix[re_total] Omega_mean_logsd = multiply_lower_tri_self_transpose(mu_logsd_betas_random_L);
}
