data {
  int<lower=1> Nd;                 // dimensions (300)
  int<lower=1> Nr;                 // traditions
  int<lower=1> Ne;                 // epigones
  int<lower=1> K;                  // number of factors

  // Epigone books
  int<lower=1> Nb;
  matrix[Nb, Nd] yE;               // each row is a book (point in Nd)
  array[Nb] int<lower=1,upper=Nr> rIDE;
  array[Nb] int<lower=1,upper=Ne> eID;

  // Pioneer books (used ONLY to learn covariance shape)
  int<lower=1> Np;
  matrix[Np, Nd] yP;
  array[Np] int<lower=1,upper=Nr> rIDP;

  // Means (computed in R from data)
  vector[Nd] muG;                  // global mean across all books
  matrix[Nr, Nd] muP;              // pioneer/tradition mean across pioneer books

  // SD predictors (computed in R from pioneer books / global)
  vector[Nd] sdd;                  // dimension-level SD predictor
  vector[Nr] sdr;                  // tradition-level RMS SD predictor
  matrix[Nr, Nd] sdp;              // tradition-by-dimension SD predictor (from pioneers)
}

parameters {
  // Attraction baseline
  real a;

  // Learning variance component (estimation error around pioneer mean)
  real<lower=0> lsb;
  real<lower=0> lnud;
  real<lower=0> lnur;
  real<lower=0> lnup;

  // Innovation variance component (style-finding)
  real<lower=0> dsb;
  real<lower=0> dnud;
  real<lower=0> dnur;
  real<lower=0> dnup;

  // Epigone within-career/book-to-book dispersion around style mean
  real<lower=0> sb;
  real<lower=0> nud;
  real<lower=0> nur;
  real<lower=0> nup;

  // Epigone-level effects: 10 params
  // 1 = avg attraction per epigone (scalar)
  // 2..4 = learning variance modifiers (for lsb, lnud, lnup)
  // 5..7 = innovation variance modifiers (for dsb, dnud, dnup)
  // 8..10 = dispersion variance modifiers (for sb, nud, nup)
  matrix[10, Ne] zV;
  vector<lower=0>[10] sigmaV;
  cholesky_factor_corr[10] L_Rho_V;

  // Dimension-specific attraction deviations around avg attraction (independent across dims here)
  matrix[Ne, Nd] zAdev;
  vector<lower=0>[Nd] sigmaAdev;

  // Learning + innovation shocks per epigone per dim (iid across dims here)
  matrix[Ne, Nd] z_le;
  matrix[Ne, Nd] z_de;

  // ===== Factor covariance (tradition-specific covariance shape) =====
  matrix[Nd, K] zLambda;
  vector<lower=0>[K] lambda_scale;

  // tradition-specific factor scales (Nr x K), centered per factor so avg scale ~ 1
  matrix[Nr, K] tau_raw;

  // latent factor scores per book
  matrix[Nb, K] fE;
  matrix[Np, K] fP;
}

transformed parameters {
  matrix[Ne, 10] V = (diag_pre_multiply(sigmaV, L_Rho_V) * zV)'; // Ne x 10
  matrix[Nd, K] Lambda = zLambda * diag_matrix(lambda_scale);

  matrix[Nr, K] tau;
  for (k in 1:K) {
    real m = mean(tau_raw[, k]);
    for (r in 1:Nr)
      tau[r, k] = exp(tau_raw[r, k] - m);   // average across r is about 1
  }
}

model {
  // Priors
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

  to_vector(zV) ~ normal(0, 1);
  sigmaV ~ exponential(1);
  L_Rho_V ~ lkj_corr_cholesky(2);

  to_vector(zAdev) ~ normal(0, 1);
  sigmaAdev ~ exponential(1);

  to_vector(z_le) ~ normal(0, 1);
  to_vector(z_de) ~ normal(0, 1);

  to_vector(zLambda) ~ normal(0, 1);
  lambda_scale ~ exponential(1);

  to_vector(tau_raw) ~ normal(0, 0.5);

  to_vector(fE) ~ normal(0, 1);
  to_vector(fP) ~ normal(0, 1);

  // =========================
  // PIONEERS: learn covariance SHAPE only
  // Standardize pioneer residuals by sdp[r,d] (computed from pioneers in R)
  // =========================
  for (p in 1:Np) {
    int r = rIDP[p];

    vector[Nd] muP_r = muP[r]';
    vector[Nd] sdP_r = sdp[r]';
    for (d in 1:Nd)
      if (sdP_r[d] < 1e-8) sdP_r[d] = 1e-8;

    vector[Nd] z = (yP[p]' - muP_r) ./ sdP_r;

    vector[K] tau_r = tau[r]';
    vector[K] f = fP[p]';
    vector[Nd] m = Lambda * (tau_r .* f);

    // residual noise is 1 on the standardized scale
    z ~ normal(m, 1);
  }

  // =========================
  // EPIGONES: mean process + dispersion scaling + same covariance shape
  // =========================
  for (b in 1:Nb) {
    int e = eID[b];
    int r = rIDE[b];

    vector[Nd] muP_r = muP[r]';
    vector[K] tau_r = tau[r]';
    vector[K] f = fE[b]';
    vector[Nd] factor_term = Lambda * (tau_r .* f);

    vector[Nd] lsig;
    vector[Nd] dsig;
    vector[Nd] sig;
    vector[Nd] ymu;

    for (d in 1:Nd) {
      // NOTE: no varying effect on the tradition-average SD term (sdr[r]) parts
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

      real attract = a + V[e,1] + sigmaAdev[d] * zAdev[e,d];
      real learn   = z_le[e,d] * lsig[d];
      real innov   = z_de[e,d] * dsig[d];

      // same algebra as your original:
      // muP + attract*(muG - muP - learn) + learn + innov
      ymu[d] = muP_r[d]
             + attract * (muG[d] - muP_r[d] - learn)
             + learn
             + innov;
    }

    // Epigone book likelihood:
    // covariance shape from factors, marginal scale from sig[d]
    // => Cov = diag(sig) * (Lambda diag(tau^2) Lambda' + I) * diag(sig)
    yE[b]' ~ normal(ymu + sig .* factor_term, sig);
  }
}

generated quantities {
  matrix[10,10] Rho_V = multiply_lower_tri_self_transpose(L_Rho_V);
}
