functions {

  // The function returns the log probability p(X | alpha) of
  // a Dirichlet-Multinomial model. We integrate out the
  // Multinomial probabilities which are not of interest in our inference.
  //
  // Parameters
  // ----------
  // vector x: Observed abundances. Sum to n
  // vector alphas: Dirichlet parameters of the form I*p_i, where I is the
  //                dispersion parameter and p_i are the mean parameters
  //
  // Returns
  // -------
  // real logpdf : The logpdf of the data
  real dirichlet_multinomial(vector x, vector alphas) {
    real n;
    real alpha0;
    real norm_factor;
    real prod_factor;
    real logpdf;

    n = sum(x);
    alpha0 = sum(alphas);

    norm_factor = lgamma(n + 1) + lgamma(alpha0) - lgamma(n + alpha0);
    prod_factor = sum(lgamma(x + alphas) - (lgamma(x + 1.0) + lgamma(alphas)));
    logpdf = norm_factor + prod_factor;

    return logpdf;

  }

} data {

  int S; // Number of species
  int N; // Number of samples
  int P; // Number of predictor variables for dispersal and drift effects
  int K; // Number of predictor variables for OTU mean effects
  int G; // Number of random effect (OTU-specific) predictor variables on dispersal and drift
  vector[N] Nt; // Abundance per sample
  matrix[N, S] abund; // OTU abundance
  matrix[N, P] X; // Design matrix for fixed dispersal effects
  matrix[N, K] W; // Design matrix for fixed OTU mean effects
  matrix[N, G] Z; // Design matrix for random OTU-specific effects on dispersal

} parameters {

  vector[P] beta; // Coefficients for dispersal and drift effects
  matrix[K, S - 1] Beta_meta; // Metacommunity OTU coefficients
  matrix[G, S - 1] Omega; // Random OTU specific effects of dispersal
  real<lower=0> sigmas[G]; // Variance for OTU-specific random effects of dispersal

} transformed parameters {

  matrix[N, S - 1] meta_p;
  matrix[P, (S - 1)] Beta; // Full fixed effects matrix

  meta_p = inv_logit(W * Beta_meta);

  // A helper matrix for easy matrix multiplication below
  for(i in 1:(S - 1)){
    Beta[:, i] = beta;
  }

} model {

  vector[2] alphas;
  vector[2] abund_OTU;
  matrix[N, S - 1] I;

  // Priors
  beta ~ normal(0, 3);
  for(g in 1:G){
    Omega[g, :] ~ normal(0, sigmas[g]); // Random effect priors
    sigmas[g] ~ normal(0, 2); // Half-normal on sigmas
  }

  // Priors on mean effects
  for(k in 1:K){
    Beta_meta[k, :] ~ normal(0, 1); // Tight constraints as some might be all 0
  }

  // Accounting for fixed and random effects on dispersal
  I = exp(X*Beta + Z*Omega);

  // Loop through OTUs
  for(i in 1:(S - 1)){

    // Loop through samples within OTUs
    for(j in 1:N){

      alphas[1] = meta_p[j, i] * I[j, i];
      alphas[2] = (1 - meta_p[j, i]) * I[j, i];
      abund_OTU[1] = abund[j, i];
      abund_OTU[2] = Nt[j] - abund[j, i];

      target += dirichlet_multinomial(abund_OTU, alphas);
    }

  }
} generated quantities{

  vector[2] alphas;
  vector[2] abund_OTU;
  matrix[N, S - 1] I;
  real log_lik[N];

  // Initialize all log_likes as zero
  for(j in 1:N)
    log_lik[j] = 0;

  I = exp(X*Beta + Z*Omega);

  // Loop through OTUs
  for(i in 1:(S - 1)){
    // Loop through samples within OTUs
    for(j in 1:N){

      alphas[1] = meta_p[j, i] * I[j, i];
      alphas[2] = (1 - meta_p[j, i]) * I[j, i];
      abund_OTU[1] = abund[j, i];
      abund_OTU[2] = Nt[j] - abund[j, i];

      log_lik[j] += dirichlet_multinomial(abund_OTU, alphas);
    }

  }
}