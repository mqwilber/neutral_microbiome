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
  vector[N] Nt; // Abundance per sample
  matrix[N, S] abund; // OTU relative abundance
  matrix[N, P] X; // Design matrix

} parameters {

  vector[P] beta; // Coefficients for the effects on dispersal and drift

  // Metacommunity relative abundances. Inferred from the data. Alternatively,
  // If these are known they could be passed in as data and inference could
  // proceed as before e.g. similar to the analysis in Burns et al. 2016 and
  // Loudon et al. 2016 where the metacommunity is assumed to be known a priori
  // This would simply requiring moving simplex[S] meta_p; into the data block
  simplex[S] meta_p;
}
model {

  vector[2] alphas; // Hold the beta-binomial shape parameters
  vector[2] abund_OTU; // Holds the abundance data for the beta-binomial
  vector[N] I; // Fundamental recruitment number (or dispersion)

  // Prior distribution om dispersal effects. Break out intercept if necessary
  beta ~ normal(0, 5);

  // Log-link effects of dispersal parameters on I
  I = exp(X*beta);

  // Loop through OTUs
  for(i in 1:(S - 1)){
    // Loop through samples within OTUs
    for(j in 1:N){

      alphas[1] = meta_p[i] * I[j];
      alphas[2] = (1 - meta_p[i]) * I[j];
      abund_OTU[1] = abund[j, i];
      abund_OTU[2] = Nt[j] - abund[j, i];

      // Note that this is a two-dimensional dirichlet-multinomial, which is
      // equivalent to a beta-binomial
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

      // Given the independent species assumption, the log-likelihood of each
      // OTU for each sample are added together.
      log_lik[j] += dirichlet_multinomial(abund_OTU, alphas);
    }

  }
}