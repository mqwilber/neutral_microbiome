functions {

  // Define a function for the collapsed dirichlet-multinomial
  // Specifically, the function returns the probability p(X | \alpha) based
  // On a Dirichlet-Multinomial parameterization. We integrate out the Multinomial
  // probabilities which we really don't care about.
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
  int P; // Number of predictor variables for dispersal effects 
  int K; // Number of predictor variables for OTU mean effects
  int G; // Number of random effect (OTU-specific) predictor variables
  vector[N] Nt; // Abundance per sample
  matrix[N, S] abund; // OTU abundance
  matrix[N, P] X; // Design matrix for fixed dispersal effects
  matrix[N, K] W; // Design matrix for fixed OTU mean effects
  matrix[N, G] Z; // Design matrix for random OTU-specific effects on dispersal

} parameters {

  vector[P] beta; // Coefficients for the link function
  matrix[K, S - 1] Beta_meta; // Metacommunity OTU coefficients
  matrix[G, S - 1] Omega; // Random OTU specific effects of dispersal
  real<lower=0> sigmas[G];

} transformed parameters {

  matrix[N, S - 1] meta_p;
  matrix[P, (S - 1)] Beta; // Full fixed effects matrix

  meta_p = inv_logit(W * Beta_meta);

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
    sigmas[g] ~ normal(0, 2);
  }

  for(k in 1:K){
    Beta_meta[k, :] ~ normal(0, 1); // Tight constraints as some might be all 0
  }

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