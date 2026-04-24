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
  
  matrix[10, Ne] zVe;
  vector<lower=0>[10] sigmVe;
  cholesky_factor_corr[10] L_Rho_e;
  
  matrix[Ne, Nd] z_le; //error variance per epigone
  matrix[Ne, Nd] z_de; //exploratory variance per epigone

}
transformed parameters {
  vector[N] mu;
  vector[N] lsig;
  vector[N] dsig;
  vector[N] sig;
  
  matrix[Ne, 10] Ve = (diag_pre_multiply(sigmVe, L_Rho_e) * zVe)';

  for (i in 1:N) {
    lsig[i] = sqrt(exp(log(lsb)+Ve[eID[i],2])^2 + exp(log(lnud)+Ve[eID[i],3])^2 * sdd[i]^2 + lnur^2 * sdr[i]^2 + exp(log(lnup)+Ve[eID[i],4])^2 * sdp[i]^2);
    dsig[i] = sqrt(exp(log(dsb)+Ve[eID[i],5])^2 + exp(log(dnud)+Ve[eID[i],6])^2 * sdd[i]^2 + dnur^2 * sdr[i]^2 + exp(log(dnup)+Ve[eID[i],7])^2 * sdp[i]^2);
    sig[i] = sqrt(exp(log(sb)+Ve[eID[i],8])^2 + exp(log(nud)+Ve[eID[i],9])^2 * sdd[i]^2 + nur^2 * sdr[i]^2 + exp(log(nup)+Ve[eID[i],10])^2 * sdp[i]^2);
    mu[i] = mup[i] + (a + Ve[eID[i],1])*(mud[i] - mup[i] - z_le[eID[i],dim[i]] * lsig[i]) + z_le[eID[i],dim[i]] * lsig[i] + z_de[eID[i],dim[i]] * dsig[i];
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
  
  to_vector(zVe) ~ normal(0, 1);
  sigmVe ~ exponential(1);
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
    matrix[10, 10] Rho_e;
    Rho_e = multiply_lower_tri_self_transpose(L_Rho_e);
}

