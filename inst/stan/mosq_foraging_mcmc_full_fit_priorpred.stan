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
  real P_FA(real alpha, real mu, real PB){
    return P_F(alpha, mu) * PB;
  }
  real P_FD(real alpha, real mu, real PB){
    return P_F(alpha, mu) * (1 - PB);
  }
}

data {
  int<lower=0> N;                // number of observations
  int<lower=0> tr;               // number of insecticides
  int treat[N];                  // treatment group index for each obs
  int<lower=0> y[N,4];           // multinomial outcomes per experiment
  int<lower=0, upper=1> use_likelihood;   // 1 = normal fit, 0 = prior predictive
}

parameters {
  vector<lower=0, upper=1>[tr+1] InitialPostprandialkillingEfficacy;
  vector<lower=0>[tr+1] KillingDuringHostSeeking;
  vector<lower=0, upper=1>[tr+1] InitialRepellencyRate;

  real a;
  real m;
  real b;
}

transformed parameters {
  real<lower=0> alpha_k0 = exp(a);
  real<lower=0> mu_k0    = exp(m);
  real<lower=0, upper=1> pB_k0 = inv_logit(b);
}

model {
  matrix[4,N] theta;

  // Priors
  target += uniform_lpdf(InitialPostprandialkillingEfficacy | 0, 1);
  target += lognormal_lpdf(1 - InitialRepellencyRate | 0, 5);
  target += lognormal_lpdf(KillingDuringHostSeeking | 0, 5);

  target += normal_lpdf(a | 0, 6);
  target += normal_lpdf(m | 0, 6);
  target += logistic_lpdf(b | 0, 1);

  // Likelihood (optional)
  for (i in 1:N) {

    if (treat[i] == 0) {
      theta[1,i] = P_UA(alpha_k0, mu_k0);
      theta[2,i] = P_UD(alpha_k0, mu_k0);
      theta[3,i] = P_FA(alpha_k0, mu_k0, pB_k0);
      theta[4,i] = P_FD(alpha_k0, mu_k0, pB_k0);

    } else {
      real alpha_i = (1 - InitialRepellencyRate[treat[i]+1]) * alpha_k0;
      real mu_i    =  mu_k0 + KillingDuringHostSeeking[treat[i]+1] * alpha_k0;
      real pBi     = pB_k0 * (1 - InitialPostprandialkillingEfficacy[treat[i]+1]);

      theta[1,i] = P_UA(alpha_i, mu_i);
      theta[2,i] = P_UD(alpha_i, mu_i);
      theta[3,i] = P_FA(alpha_i, mu_i, pBi);
      theta[4,i] = P_FD(alpha_i, mu_i, pBi);
    }

    if (use_likelihood == 1)
      target += multinomial_lpmf(y[i,] | theta[,i]);
  }
}

generated quantities {
  vector[tr+1] InitialPreprandialkillingEfficacy;

  real<lower=0> alpha_0 = exp(a);
  real<lower=0> mu_0    = exp(m);
  real<lower=0, upper=1> pB_0 = inv_logit(b);

  vector<lower=0>[tr+1] alpha_i = (1 - InitialRepellencyRate) * alpha_0;
  vector<lower=0>[tr+1] mu_i    = mu_0 + KillingDuringHostSeeking * alpha_0;

  real<lower=0, upper=1> PSc =
    1 - ((1 - exp(-(alpha_0 + mu_0))) * mu_0 / (alpha_0 + mu_0));

  vector<lower=0, upper=1>[tr+1] PS;

  InitialPreprandialkillingEfficacy[1] = 0;
  PS[1] = PSc;

  for (j in 1:tr) {
    PS[j+1] =
      1 - ((1 - exp(-(alpha_i[j+1] + mu_i[j+1]))) *
            mu_i[j+1] / (alpha_i[j+1] + mu_i[j+1]));

    InitialPreprandialkillingEfficacy[j+1] =
      1 - (PS[j+1] / PSc);
  }
}
