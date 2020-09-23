functions {
  int[] compute_n_k(int K, int[] group) {
    int n_k[K] = rep_array(0, K);
    for(i in 1:size(group)) {
      n_k[group[i]] += 1;
    }
    return(n_k);
  }

  int seq_from_to(int a, int b) {
    int n = b - a + 1;
    int out[n];
    for(i in 1:n) {
      out[i] = i;
    }
    out += a - 1;
    return(out);
  }
}

data {
  int N; // (Total) Observations
  int K; // Groups
  int F; // Outcome variables
  int P; // Location predictors
  int Q; // Scale predictors
  // Stay simple here; no randoms.
  /* int P_random; // Location RE slopes */
  /* int Q_random; // Scale RE slopes */
  /* int P_random_ind[P_random]; // Indices in x_loc for RE location slopes */
  /* int Q_random_ind[Q_random]; // Indices in x_sca for RE scale slopes */
  int AR_P; // Lags for AR location component
  int MA_P; // Lags for MA location component
  int AR_Q; // Lags for AR scale component (?)
  int MA_Q; // Lags for MA scale component (?)

  int group[N];
  matrix[N, P] x_loc;
  matrix[N, Q] x_sca;
  matrix[N, F] y;

  /* If sorted, this is not necessary. Otherwise, would need it, and need to create an index matrix. */
  /* int time[N]; // Time point data; per observation. ASSUME y is sorted by time within K. */

  // Options
  int<lower=0, upper=1> prior_only;
}

transformed data {
  int n_k[K] = compute_n_k(K, group); // Number of time points per group (ARMA).
}

parameters {
  row_vector[F] nu; // Intercepts
  row_vector[F] sigma; // Logsd-intercepts
  matrix[P, F] mu_beta; // Location coefficients
  matrix[Q, F] logsd_beta; // Scale coefficients

  // ARMA parameters for location and scale
  // For now, no cross-lags.
  // When REs are introduced: This is going to be hard. Each group's AR coef must be [0,1]; 
  // May want to use a logistic-RE model then; ar_coef[k] = inv_logit(fe + u[k])
  // And when cross-lags are allowed, then the SUM must be [0, 1], I think.
  // Unsure about constraint on the log-sd models though; those are multiplicative.
  /* matrix[F, F] ar_location[AR_P]; */
  row_vector<lower=0, upper=1>[F] ar_location[AR_P];
  row_vector<lower=0, upper=1>[F] ma_location[MA_P];
  row_vector[F] ar_scale[AR_Q];
  row_vector[F] ma_scale[MA_Q];
  matrix[N, F] var_innovation_z;
  row_vector<lower=0>[F] var_innovation_sd;
  /* cholesky_factor_corr[F] var_innovation_cor; */ // Test this later!

  // Error Cor
  cholesky_factor_corr[F] epsilon_cor_L;

  // Auxilliary parameters
  
}

transformed parameters {
  matrix[N, F] var_innovation = diag_post_multiply(var_innovation_z, var_innovation_sd);

  // FE
  matrix[N, F] mu_hat = rep_matrix(nu, N) + x_loc * mu_beta;
  matrix[N, F] logsd_hat = rep_matrix(sigma, N) + x_sca * logsd_beta;
  logsd_hat += var_innovation;
  // RE

  // Innovations

  // ARMA Location (No cross-lags, so element-wise coef multiplication)
  for(n in 2:N) {
    // AR
    for(lag in 1:min(AR_P, n - 1)) {
      mu_hat[n] += y[n-lag] .* ar_location[lag];
    }
    // MA
    for(lag in 1:min(MA_P, n - 1)) {
      mu_hat[n] += (y[n - lag] - mu_hat[n - lag]) .* ma_location[lag];
    }
  }
  // ARMA Scale
  for(n in 2:N) {
    // AR
    for(lag in 1:min(AR_Q, n - 1)) {
      logsd_hat[n] += logsd_hat[n - lag] .* ar_scale[lag];
    }
    // MA
    for(lag in 1:min(MA_Q, n - 1)) {
      logsd_hat[n] += var_innovation[n - lag] .* ma_scale[lag];
    }
  }
}

model {
  nu ~ std_normal();
  sigma ~ std_normal();
  mu_beta ~ std_normal();
  logsd_beta ~ std_normal();

  ar_location ~ std_normal();
  ma_location ~ std_normal();
  ar_scale ~ std_normal();
  ma_scale ~ std_normal();

  var_innovation_z ~ std_normal();
  var_innovation_sd ~ std_normal();
  epsilon_cor_L ~ lkj_corr_cholesky(1);

  for(n in 1:N) {
    y[n] ~ multi_normal_cholesky(mu_hat[n], diag_pre_multiply(exp(logsd_hat[n]), epsilon_cor_L));
  }
}
