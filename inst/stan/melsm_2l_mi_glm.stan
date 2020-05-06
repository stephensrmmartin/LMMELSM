/* 
   MELSM-Latent model:
   X[i] = Nu + Eta[i] * Lambda + Error[i];
Eta[i] ~ N(alpha[k_i], Sigma_eta[k_i]);
alpha[k_i] ~ N(0, Cor(alpha));
Sigma_eta[k_i] = L * SD[k_i] * SD[k_i]' * L'
log(SD[k_i]) ~ N(0, ??)
 */ 
functions {
  matrix lambda_mat(int J, int F, int[] J_f, int[,] F_ind, vector lambda_est) {
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
  
}

data {
  // Meta-data
  int N; // Total number of observations.
  int J; // Total number of indicators.
  int F; // Number of latent factors
  int K; // Number of grouping variables

  int group[N]; // Grouping indicator.

  // Indicators
  int J_f[F]; // Number of indicators for each factor.
  int F_ind[F,J]; // Indicator indices.
  matrix[N, J] x;

}

transformed data {
  int lambda_total = sum(J_f);
}

parameters {
  // Measurement model
  row_vector[J] nu;
  vector<lower=0>[lambda_total] lambda_est;
  row_vector<lower=0>[J] sigma;

  // Latent vars
  matrix[N, F] eta_random_z; // eta[i] = alpha[k_i] + eta_random_z[i]*(Diag(logsd[k_i])*L_eta_random)'
  matrix[K, F*2] eta_mean_logsd_z; // [loc, logsd] = [0, 0] + eta_mean_logsd_z*(Diag(sigma_mean_logsd)*L_mean_logsd)'
  cholesky_factor_corr[F] L_eta_random; // Assumes eta cors are same across K.
  cholesky_factor_corr[F*2] L_mean_logsd;
  vector<lower=0>[F*2] sigma_mean_logsd; // No scale model on L2 SDs.
  
}


transformed parameters {
  matrix[F, J] lambda = lambda_mat(J, F, J_f, F_ind, lambda_est);
  matrix[K, F*2] eta_mean_logsd = eta_mean_logsd_z * (diag_pre_multiply(sigma_mean_logsd, L_mean_logsd))'; // both assumed to have means of zero; implicitly think [0,0,...,0] + ...
  matrix[K, F] eta_mean = eta_mean_logsd[,1:F];
  matrix[K, F] eta_sd = exp(eta_mean_logsd[,(F+1) : (F*2)]);
  matrix[N, F] eta = eta_mean[group,];
  matrix[F,F] U_group_eta[K]; // Will be Upper (transposed) Cholesky covariance
  for(k in 1:K){
    U_group_eta[k] = (diag_pre_multiply(eta_sd[k], L_eta_random))';
  }
  for(n in 1:N){
    eta[n] += eta_random_z[n] * U_group_eta[group[n]];
  }


  
}

model {

  // Priors
  nu ~ std_normal();
  lambda_est ~ std_normal();
  sigma ~ std_normal();

  to_vector(eta_random_z) ~ std_normal();
  to_vector(eta_mean_logsd_z) ~ std_normal();
  L_eta_random ~ lkj_corr_cholesky(1);
  L_mean_logsd ~ lkj_corr_cholesky(1);
  sigma_mean_logsd ~ std_normal();

  // to_vector(x) ~ normal(to_vector(xhat), to_vector(shat));
  for(j in 1:J) {
    x[,j] ~ normal_id_glm(eta, nu[j], lambda[,j], sigma[j]);
  }

}

generated quantities {
  corr_matrix[F] Omega_eta = multiply_lower_tri_self_transpose(L_eta_random);
  corr_matrix[F*2] Omega_mean_logsd = multiply_lower_tri_self_transpose(L_mean_logsd);
  
}

