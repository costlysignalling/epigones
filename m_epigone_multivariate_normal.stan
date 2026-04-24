data {
  int<lower=1> Nb;
  int<lower=1> Nd;
  int<lower=1> Nr;
  int<lower=1> Ne;

  array[Nb] int<lower=1,upper=Nr> rID;
  array[Nb] int<lower=1,upper=Ne> eID;

  matrix[Nb, Nd] y;          // each row is a book (300D)

  vector[Nd] muG;            // global mean
  matrix[Nr, Nd] muP;        // pioneer mean per tradition

  vector[Nd] sdd;            // dimension SD predictor
  vector[Nr] sdr;            // pioneer-average SD predictor (RMS)
  matrix[Nr, Nd] sdp;        // pioneer-by-dim SD predictor
}

parameters {
  real a;

  // learning variance component
  real<lower=0> lsb;
  real<lower=0> lnud;
  real<lower=0> lnur;
  real<lower=0> lnup;

  // innovation variance component
  real<lower=0> dsb;
  real<lower=0> dnud;
  real<lower=0> dnur;
  real<lower=0> dnup;

  // within-book variance component
  real<lower=0> sb;
  real<lower=0> nud;
  real<lower=0> nur;
  real<lower=0> nup;

  // dimension-specific attraction deviations (300D), correlated across dims
  matrix[Nd, Ne] zA;
  vector<lower=0>[Nd] sigmaA;
  cholesky_factor_corr[Nd] L_Rho_A;

  // 10 epigone-level random effects: avg attraction + 9 variance modifiers
  matrix[10, Ne] zV;
  vector<lower=0>[10] sigmaV;
  cholesky_factor_corr[10] L_Rho_V;

  // learning + innovation shocks per epigone per dim (still iid here)
  matrix[Ne, Nd] z_le;
  matrix[Ne, Nd] z_de;

  // within-book correlation across dims (true MVN likelihood)
  cholesky_factor_corr[Nd] L_Rho_y;
}

transformed parameters {
  matrix[Ne, Nd] Adev = (diag_pre_multiply(sigmaA, L_Rho_A) * zA)';  // Ne x Nd
  matrix[Ne, 10] V    = (diag_pre_multiply(sigmaV, L_Rho_V) * zV)';  // Ne x 10
}

model {
  // priors (keep/adjust as you like)
  a ~ normal(0, 1);

  lsb  ~ exponential(1);
  lnud ~ normal(0.5, 1);
  lnur ~ normal(0.5, 1);
  lnup ~ normal(0.5, 1);

  dsb  ~ exponential(1);
  dnud ~ normal(0.5, 1);
  dnur ~ normal(0.5, 1);
  dnup ~ normal(0.5, 1);

  sb   ~ exponential(1);
  nud  ~ normal(0.5, 1);
  nur  ~ normal(0.5, 1);
  nup  ~ normal(0.5, 1);

  to_vector(zA) ~ normal(0, 1);
  sigmaA ~ exponential(1);
  L_Rho_A ~ lkj_corr_cholesky(8);    // strong shrinkage in 300D

  to_vector(zV) ~ normal(0, 1);
  sigmaV ~ exponential(1);
  L_Rho_V ~ lkj_corr_cholesky(2);

  to_vector(z_le) ~ normal(0, 1);
  to_vector(z_de) ~ normal(0, 1);

  L_Rho_y ~ lkj_corr_cholesky(10);   // strong shrinkage for book covariance

  // likelihood
  for (b in 1:Nb) {
    int e = eID[b];
    int r = rID[b];

    vector[Nd] lsig;
    vector[Nd] dsig;
    vector[Nd] sig;
    vector[Nd] muP_r = to_vector(muP[r]');
    vector[Nd] ymu;

    for (d in 1:Nd) {
      // corrected varying-effect allocation:
      // - sdr[r] terms have NO V[e,*]
      real vL = exp(log(lsb)  + V[e,2])
              + exp(log(lnud) + V[e,3]) * square(sdd[d])
              + exp(log(lnur))          * square(sdr[r])
              + exp(log(lnup) + V[e,4]) * square(sdp[r,d]);

      real vD = exp(log(dsb)  + V[e,5])
              + exp(log(dnud) + V[e,6]) * square(sdd[d])
              + exp(log(dnur))          * square(sdr[r])
              + exp(log(dnup) + V[e,7]) * square(sdp[r,d]);

      real vB = exp(log(sb)   + V[e,8])
              + exp(log(nud)  + V[e,9]) * square(sdd[d])
              + exp(log(nur))           * square(sdr[r])
              + exp(log(nup)  + V[e,10]) * square(sdp[r,d]);

      lsig[d] = sqrt(vL);
      dsig[d] = sqrt(vD);
      sig[d]  = sqrt(vB);

      // attraction: global + epigone-average + epigone-dim deviation
      real attract = a + V[e,1] + Adev[e,d];

      real learn = z_le[e,d] * lsig[d];
      real innov = z_de[e,d] * dsig[d];

      ymu[d] = muP_r[d]
             + attract * (muG[d] - muP_r[d] - learn)
             + learn
             + innov;
    }

    matrix[Nd, Nd] L_Sigma = diag_pre_multiply(sig, L_Rho_y);
    y[b]' ~ multi_normal_cholesky(ymu, L_Sigma);
  }
}

generated quantities {
  matrix[Nd, Nd] Rho_y = multiply_lower_tri_self_transpose(L_Rho_y);
  matrix[Nd, Nd] Rho_A = multiply_lower_tri_self_transpose(L_Rho_A);
  matrix[10,10]  Rho_V = multiply_lower_tri_self_transpose(L_Rho_V);
}
