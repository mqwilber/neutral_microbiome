## Stan models for fitting flexible community assembly models

- `neutral_beta_binomial.stan`: Neutral community assembly model with beta-binomial sampling to account for sampling processes (see Appendix S1 of the manuscript).  This model assumes a constant metacommunity with no selection processes and allows dispersal and drift processes to vary among local communities and the community scale.  This is model can be used to test Hypothesis 1: Base neutral model or Hypothesis 2: Dispersal and drift vary across local communities.
    - The data this model needs are described in the Stan model. In short
        1. `S`: The number of species that are being considered. Only S - 1 species are fit. The Sth species is ``other`` consisting of all other species.
        1. `N`: The number of local community samples. Each sample should be associated with a vector of community relative abundances and any relevant covariates.
        1. `Nt`: Sampled abundance for each `N` samples. This is used because all NT individuals in a community can rarely be sampled for microbial communities and Nt is used is the dirichlet-multinomial sampling model.
        1. `abund`: The N (relative) abundance vectors of each community.
        1. `X`: A design matrix that holds the covariates that may affect community-level dispersal and drift.
        1. `P`: The number of predictor variables in `X`.
- `neutral_beta_binomial_meaneffects.stan`: A community assembly model with beta-binomial sampling to account for community-level sampling processes (see Appendix S1 of the manuscript). This model is similar to `neutral_beta_binomial.stan`, but allows for species-specific selection processes in addition to community-level effects of dispersal and drift.  This model can be used to test variants of Hypothesis 1 and Hypothesis 2, where the metacommunity is allowed to vary by location.
    - This model has the same parameters as the previous model with two additions.
        1. `W`: A design matrix that includes covariates affecting species-specific selection effects.
        1. `K`: The number of predictor variables in `W`.
- `neutral_beta_binomial_meanrandomeffects.stan`: A community assembly model with beta-binomial sampling. This model includes species-specific random effects of dispersal and drift and can be used to test Hypothesis 3: Dispersal and drift vary by bacterial species and Hypothesis 4: Selection processes drive community assembly.
    - This model has the same parameters as the above model with the additional parameters
        1. `Z`: A design matrix for species-specific effects on dispersal and drift.
        1. `G`: The number of predictor variables in `Z`.
- `neutral_dirichlet_multinomial_meaneffects.stan`: This is the same model as `neutral_beta_binomial_meaneffects.stan`, but does not make an independent species assumption. Because of this, this model cannot be used to test Hypothesis 3, but can be used to compare whether inference on selection and community-level dispersal and drift effects changes under the independent species assumption.




