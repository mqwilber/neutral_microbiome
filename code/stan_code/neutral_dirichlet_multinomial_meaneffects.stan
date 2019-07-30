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
  int P; // Number of predictor variables for dispersal effects
  int K; // Number of predictor variables for OTU mean effects
  real abund[N, S]; // OTU abundance
  matrix[N, P] X; // Design matrix for fixed dispersal effects
  matrix[N, K] W; // Design matrix for fixed OTU mean effects

} parameters {

  vector[P] beta; // Fixed effect dispersal coefficients
  matrix[K, S - 1] Beta_meta; // Metacommunity OTU coefficients

} transformed parameters {

  simplex[S] meta_p[N]; // Array of size N with simplexes with S elements
  matrix[N, S - 1] meta_p_logit;
  real sum_meta;

  meta_p_logit = W * Beta_meta;

  // Compute meta_pusing a softmax link. For example, see page 97 in
  // Extending the Linear Model with R, Faraway 2007.
  for(i in 1:N){
    sum_meta = 1 + sum(exp(meta_p_logit[i, :]));
    for(s in 1:S) {
        if(s != S)
          meta_p[i][s] = exp(meta_p_logit[i, s]) / sum_meta;
        else
          meta_p[i][s] = 1 - sum(meta_p[i][1:(S - 1)]);
    }
  }

} model {

  vector[S] alphas; // Dirichlet parameters of the neutral model
  vector[N] I; // Fundamental recruitment number

  // Priors on the dispersal parameters
  for(p in 1:P) {
    if(p == 1)
      beta[p] ~ normal(0, 3); // Adjust the prior on the intercept if desired
    else
      beta[p] ~ normal(0, 3);
  }

  // Priors on the OTU-specific mean effects
  for(k in 1:K){
    if(k == 1)
      Beta_meta[k, :] ~ normal(0, 5); // intercepts
    else
      Beta_meta[k, :] ~ normal(0, 2);
  }

  // Varying dispersion by sample. Modeled on the log-link scale
  I = exp(X*beta);

  // Loop through samples. Each sample is a vector of observed abundances
  for(j in 1:N){

    alphas = I[j] * meta_p[j]; // Dispersion * mean
    target += dirichlet_multinomial(to_vector(abund[j]), to_vector(alphas));

  }

} generated quantities {

  vector[S] alphas;
  vector[N] I;
  vector[N] log_lik;

  I = exp(X*beta);

  // Loop through samples within OTUs
  for(j in 1:N) {

    alphas = I[j] * meta_p[j]; // Dispersion * mean
    log_lik[j] = dirichlet_multinomial(to_vector(abund[j]), to_vector(alphas));

  }

}
