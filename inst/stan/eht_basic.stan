// Champagne et al. (2025) Nature Communications 16:11162, https://doi.org/10.1038/s41467-025-66130-y
// Stan model for experimental hut trials (EHT), control and intervention arms
// Outcomes per observation: unfed alive (UA), unfed dead (UD), fed alive (FA), fed dead (FD)
// pC is the probability of surviving after feeding, i.e. the product P_C P_D P_E of
// Chitnis et al. (2008) (finding a resting place, resting, ovipositing), not P_B (biting)

functions{
  // probabilities in the Markov chain model for mosquito foraging behaviour
  real P_UA(real alpha, real mu){
    return exp(-(alpha + mu));
  }
  real P_F(real alpha, real mu){
    return (1 - P_UA(alpha, mu)) * alpha / (alpha + mu);
  }
  real P_UD(real alpha, real mu){
    return (1 - P_UA(alpha, mu)) * mu / (alpha + mu);
  }
  real P_FA(real alpha, real mu, real PC){
    return P_F(alpha, mu) * PC;
  }
  real P_FD(real alpha, real mu, real PC){
    return P_F(alpha, mu) * (1 - PC);
  }
}

data {
  int<lower=0> N;                // number of observations
  int<lower=0> tr;               // number of insecticides
  array[N] int treat;            // treatment group index for each obs
  array[N, 4] int<lower=0> y;    // multinomial outcomes per experiment
  int<lower=0, upper=1> use_likelihood;   // 1 = normal fit, 0 = prior predictive
}

parameters {
  vector<lower=0, upper=1>[tr+1] PostprandialKilling;
  vector<lower=0>[tr+1] KillingDuringHostSeeking;
  vector<lower=0, upper=1>[tr+1] HostAvailabilityReduction;

  real a;
  real m;
  real b;
}

transformed parameters {
  real<lower=0> alpha_k0 = exp(a);
  real<lower=0> mu_k0    = exp(m);
  real<lower=0, upper=1> pC_k0 = inv_logit(b);
}

model {
  matrix[4,N] theta;

  // Priors
  target += uniform_lpdf(PostprandialKilling | 0, 1);
  target += lognormal_lpdf(1 - HostAvailabilityReduction | 0, 5);
  target += lognormal_lpdf(KillingDuringHostSeeking | 0, 5);

  target += normal_lpdf(a | 0, 6);
  target += normal_lpdf(m | 0, 6);
  target += logistic_lpdf(b | 0, 1);

  // Likelihood (optional)
  for (i in 1:N) {

    if (treat[i] == 0) {
      theta[1,i] = P_UA(alpha_k0, mu_k0);
      theta[2,i] = P_UD(alpha_k0, mu_k0);
      theta[3,i] = P_FA(alpha_k0, mu_k0, pC_k0);
      theta[4,i] = P_FD(alpha_k0, mu_k0, pC_k0);

    } else {
      real alpha_i = (1 - HostAvailabilityReduction[treat[i]+1]) * alpha_k0;
      real mu_i    =  mu_k0 + KillingDuringHostSeeking[treat[i]+1] * alpha_k0;
      real pCi     = pC_k0 * (1 - PostprandialKilling[treat[i]+1]);

      theta[1,i] = P_UA(alpha_i, mu_i);
      theta[2,i] = P_UD(alpha_i, mu_i);
      theta[3,i] = P_FA(alpha_i, mu_i, pCi);
      theta[4,i] = P_FD(alpha_i, mu_i, pCi);
    }

    if (use_likelihood == 1)
      target += multinomial_lpmf(y[i,] | theta[,i]);
  }
}

generated quantities {
  vector[tr+1] PreprandialKilling;

  real<lower=0> alpha_0 = exp(a);
  real<lower=0> mu_0    = exp(m);
  real<lower=0, upper=1> pC_0 = inv_logit(b);

  vector<lower=0>[tr+1] alpha_i = (1 - HostAvailabilityReduction) * alpha_0;
  vector<lower=0>[tr+1] mu_i    = mu_0 + KillingDuringHostSeeking * alpha_0;

  real<lower=0, upper=1> PSc =
    1 - ((1 - exp(-(alpha_0 + mu_0))) * mu_0 / (alpha_0 + mu_0));

  vector<lower=0, upper=1>[tr+1] PS;

  // replicated data for posterior predictive checks
  array[N, 4] int y_rep;

  PreprandialKilling[1] = 0;
  PS[1] = PSc;

  for (j in 1:tr) {
    PS[j+1] =
      1 - ((1 - exp(-(alpha_i[j+1] + mu_i[j+1]))) *
            mu_i[j+1] / (alpha_i[j+1] + mu_i[j+1]));

    PreprandialKilling[j+1] =
      1 - (PS[j+1] / PSc);
  }

  // same probabilities as in the model block, same number of mosquitoes per observation
  for (i in 1:N) {
    vector[4] theta_rep;
    if (treat[i] == 0) {
      theta_rep = [P_UA(alpha_0, mu_0), P_UD(alpha_0, mu_0),
                   P_FA(alpha_0, mu_0, pC_0), P_FD(alpha_0, mu_0, pC_0)]';
    } else {
      real pCi = pC_0 * (1 - PostprandialKilling[treat[i]+1]);
      theta_rep = [P_UA(alpha_i[treat[i]+1], mu_i[treat[i]+1]),
                   P_UD(alpha_i[treat[i]+1], mu_i[treat[i]+1]),
                   P_FA(alpha_i[treat[i]+1], mu_i[treat[i]+1], pCi),
                   P_FD(alpha_i[treat[i]+1], mu_i[treat[i]+1], pCi)]';
    }
    y_rep[i] = multinomial_rng(theta_rep, sum(y[i]));
  }
}
