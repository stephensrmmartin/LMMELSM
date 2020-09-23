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
  /* matrix[F, F] ar_location[AR_P]; */
  row_vector[F] ar_location[AR_P];
  row_vector[F] ma_location[MA_P];
  row_vector[F] ar_scale[AR_Q];
  row_vector[F] ma_scale[MA_Q];

  // Error Cor
  cholesky_factor_corr[F] epsilon_cor;

  // Auxilliary parameters
  
}

