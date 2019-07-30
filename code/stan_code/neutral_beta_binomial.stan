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
  int P; // Predictor variables
  vector[N] Nt; // Abundance per sample
  matrix[N, S] abund; // OTU relative abundance
  matrix[N, P] X; // Design matrix
    
} parameters {

  vector[P] beta; // Coefficients for the link function
  simplex[S] meta_p; // Metacommunity relative abundances 

}  
model {

  vector[2] alphas;
  vector[2] abund_OTU;
  vector[N] I;

  beta ~ normal(0, 5);

  I = exp(X*beta);
  
  // Loop through OTUs
  for(i in 1:(S - 1)){
    // Loop through samples within OTUs
    for(j in 1:N){ 

      alphas[1] = meta_p[i] * I[j];
      alphas[2] = (1 - meta_p[i]) * I[j];
      abund_OTU[1] = abund[j, i];
      abund_OTU[2] = Nt[j] - abund[j, i];

      target += dirichlet_multinomial(abund_OTU, alphas);
    }

  }
} generated quantities{

  vector[2] alphas;
  vector[2] abund_OTU;
  vector[N] I;
  real log_lik[N];

  // Initialize all log_likes as zero
  for(j in 1:N)
    log_lik[j] = 0;

  I = exp(X*beta);

  // Loop through OTUs
  for(i in 1:(S - 1)){
    // Loop through samples within OTUs
    for(j in 1:N){ 

      alphas[1] = meta_p[i] * I[j];
      alphas[2] = (1 - meta_p[i]) * I[j];
      abund_OTU[1] = abund[j, i];
      abund_OTU[2] = Nt[j] - abund[j, i];

      log_lik[j] += dirichlet_multinomial(abund_OTU, alphas);
    }

  }
}