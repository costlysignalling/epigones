data {
  int<lower=1> N;
  int<lower=1> Nd;
  int<lower=1> Nr;
  int<lower=1> Ne;
  array[N] int rID;
  array[N] int eID;
  array[N] int book;
  array[N] int dim;
  vector[N] epi;
  vector[N] mud;
  vector[N] mup;
  vector[N] sdp;
  vector[N] sdr;
  vector[N] sdd;
}
parameters {
  real a;

  real<lower=0> lsb;
  real<lower=0> lnud;
  real<lower=0> lnur;
  real<lower=0> lnup;
  
  real<lower=0> dsb;
  real<lower=0> dnud;
  real<lower=0> dnur;
  real<lower=0> dnup;
  
  real<lower=0> sb;
  real<lower=0> nud;
  real<lower=0> nur;
  real<lower=0> nup;
  
  matrix[Nd, Ne] za_e;
  vector<lower=0>[Nd] sigma_e;
  cholesky_factor_corr[Nd] L_Rho_e;
  
  matrix[Ne, Nd] z_le; //error variance per epigone
  matrix[Ne, Nd] z_de; //exploratory variance per epigone

}
transformed parameters {
  vector[N] mu;

  vector[N] lsig = sqrt(lsb^2 + lnud^2 * sdd^2 + lnur^2 * sdr^2 + lnup^2 * sdp^2);
  vector[N] dsig = sqrt(dsb^2 + dnud^2 * sdd^2 + dnur^2 * sdr^2 + dnup^2 * sdp^2);
  vector[N] sig  = sqrt(sb^2  + nud^2 * sdd^2 + nur^2 * sdr^2 + nup^2 * sdp^2);
  
  matrix[Ne, Nd] a_e = (diag_pre_multiply(sigma_e, L_Rho_e) * za_e)';

  for (i in 1:N) {
    mu[i] = mup[i] + (a + a_e[eID[i],dim[i]])*(mud[i] - mup[i] - z_le[eID[i],dim[i]] * lsig[i]) + z_le[eID[i],dim[i]] * lsig[i] + z_de[eID[i],dim[i]] * dsig[i];
  }
}
model {
  a ~ normal(0,1);
  
  lsb ~ exponential(1);
  lnud ~ normal(0.5, 1);
  lnur ~ normal(0.5, 1);
  lnup ~ normal(0.5, 1);
  
  dsb ~ exponential(1);
  dnud ~ normal(0.5, 1);
  dnur ~ normal(0.5, 1);
  dnup ~ normal(0.5, 1);
  
  sb ~ exponential(1);
  nud ~ normal(0.5, 1);
  nur ~ normal(0.5, 1);
  nup ~ normal(0.5, 1);
  
  to_vector(za_e) ~ normal(0, 1);
  sigma_e ~ exponential(1);
  L_Rho_e ~ lkj_corr_cholesky(2);

  to_vector(z_le) ~ normal(0, 1);
  to_vector(z_de) ~ normal(0, 1);
  
  // Likelihood book
  epi ~ normal(mu, sig);
}
generated quantities {
    vector[N] log_lik;
    for (i in 1:N) {
        log_lik[i] = normal_lpdf(epi[i] | mu[i], sig[i]);
    }
    matrix[Nd, Nd] Rho_e;
    Rho_e = multiply_lower_tri_self_transpose(L_Rho_e);
}

