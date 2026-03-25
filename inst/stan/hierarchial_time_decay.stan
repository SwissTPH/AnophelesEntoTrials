// from https://git.scicore.unibas.ch/idm/countrymodelling/irsparameterisation/-/tree/main/inst/stan?ref_type=heads Denz_model.stan

functions{
  // probabilities in the Markov chain model for mosquito foraging behaviour
  real P_A( real alpha, real mu){
    return exp( -(alpha+mu)  );
  }
  real P_F(real alpha, real mu){
    return (1-P_A(alpha, mu)) *alpha /(alpha +mu);
  }
  real P_M(real alpha, real mu){
    return (1-P_A(alpha, mu))*mu /(alpha +mu);
  }

}


// The input data is a vector 'y' of length 'N'-
data {
  int<lower=0> N; //number of observations
  int<lower=0> tr; //number of different insecticide
  int<lower=0> d; //number of days
  int total[N];
  int datej[N];
  int day[N];
  int treat[N];
  int fed[N];
  int FA[N];
  int<lower=0> y[N,3]; // observed outcome vectors (UA,F, UD) for n experiments from control arm
  real<lower=0> priorsigma_mean_logrates; // sigma of the normal priors on a and m
  real<lower=0> hierarchy; // variance for phi and psi, limit --> 0 means eliminating hierarchy
}

// The parameters accepted by the model-
parameters{
   // our priors per insecticide
  vector<lower=0, upper=1>[tr+1] InitialPostprandialkillingEfficacy;
  vector<lower=0>[tr+1] KillingDuringHostSeeking;
  vector<lower=0, upper=1>[tr+1] InitialRepellencyRate;
  vector<lower=0>[tr+1] beta;
  vector<lower=0>[tr+1] kappa;

  real a;
  real m;
  real r_prob_surviving_feeding_d;

  real<lower=0> sigma_a;
  real<lower=0> sigma_m;
  real<lower=0> sigma_prob_surviving_feeding_d;

  vector[d] phi;
  vector[d] psi;
  vector[d] omega_prob_surviving_feeding_d;
}


transformed parameters{
  // declare: control rates
  vector<lower=0, upper=1>[d] control_prob_surviving_feeding;
  vector<lower=0, upper=1>[d] alpha_H_k0;
  vector<lower=0>[d] mu_k0;


  // define: control rates
  alpha_H_k0 = exp(a + phi *sigma_a);
  mu_k0 = exp(m + psi *sigma_m);
  control_prob_surviving_feeding = inv_logit(r_prob_surviving_feeding_d + omega_prob_surviving_feeding_d * sigma_prob_surviving_feeding_d);
}


model{
//declare local variables
  matrix[3,N] theta;
  vector[N] pB;
  real decay;

// Priors
  target += uniform_lpdf(InitialPostprandialkillingEfficacy | 0, 1);
  target += uniform_lpdf(KillingDuringHostSeeking | 0, 1);
  target += uniform_lpdf(InitialRepellencyRate | 0, 1);
 // target += normal_lpdf(log_RA| 0, 0.001);
 target += normal_lpdf(beta| 0, 1);
 target += normal_lpdf(kappa| 1, 1);

//priors
  // priors on means of logrates
  a ~ normal(0,priorsigma_mean_logrates);
  m ~ normal(0,priorsigma_mean_logrates);
  target += logistic_lpdf(r_prob_surviving_feeding_d |0,1);
  // priors on phi and psi (normalised deviation from mean of logrates for each group)
  phi ~ normal(0,1);
  psi ~ normal(0,1);
   target += normal_lpdf( omega_prob_surviving_feeding_d  | 0, 1);
  // priors on sigma_a and sigma_m (variation of deviation from mean of logrates for each group)
  target += 2*cauchy_lpdf(sigma_a | 0,hierarchy);
  target += 2*cauchy_lpdf(sigma_m | 0,hierarchy);
  target += 2*cauchy_lpdf(sigma_prob_surviving_feeding_d | 0,1);
 // loop for defining model

  for (i in 1:N) {

    if (treat[i]==0) {
    pB[i]   = control_prob_surviving_feeding[day[i]];
    #' removed for all theta
    theta[1,i] = (P_A(alpha_H_k0[day[i]], mu_k0[day[i]]));
    theta[2,i] = (P_F(alpha_H_k0[day[i]], mu_k0[day[i]]));
    theta[3,i] = (P_M(alpha_H_k0[day[i]], mu_k0[day[i]]));

    }
    else{
   decay= exp(-(beta[treat[i]+1]*datej[i])^kappa[treat[i]+1]);
    pB[i]   = control_prob_surviving_feeding[day[i]]*(1-InitialPostprandialkillingEfficacy[treat[i]+1]*decay);
    #' removed for all theta
    theta[1,i] = (P_A((1 -InitialRepellencyRate[treat[i]+1]*decay)*alpha_H_k0[day[i]], mu_k0[day[i]] + KillingDuringHostSeeking[treat[i]+1]*decay * alpha_H_k0[day[i]]));
    theta[2,i] = (P_F((1 -InitialRepellencyRate[treat[i]+1]*decay)*alpha_H_k0[day[i]],  mu_k0[day[i]] + KillingDuringHostSeeking[treat[i]+1]*decay * alpha_H_k0[day[i]]));
    theta[3,i] = (P_M((1 -InitialRepellencyRate[treat[i]+1]*decay)*alpha_H_k0[day[i]],  mu_k0[day[i]] + KillingDuringHostSeeking[treat[i]+1]*decay * alpha_H_k0[day[i]]));
    }
  target += binomial_lpmf(FA[i] | fed[i], pB[i]);
  target += multinomial_lpmf( y[i,] |theta[,i]);
  }


}

generated quantities {
  // for model selection with loo


vector<lower=0>[tr+1] L;
vector[tr+1] InitialPreprandialkillingEfficacy;
vector[tr+1] InitialRepellentEfficacy;

// declare: means of rates
real<lower=0> alpha_0; //
real<lower=0> mu_0; //
vector<lower=0>[tr+1] alpha_i; //
vector<lower=0>[tr+1] mu_i; //
real<lower=0, upper=1> Pfc; //
real<lower=0, upper=1> PAttc; //
vector<lower=0, upper=1>[tr+1] Pf; //
vector<lower=0, upper=1>[tr+1] PAtt; //

//transform: means of rates
alpha_0 = exp(a + sigma_a^2/2);
mu_0 = exp(m + sigma_m^2/2);
alpha_i = (1 -InitialRepellencyRate) *alpha_0;
mu_i = mu_0 + KillingDuringHostSeeking * alpha_0;
Pfc=(1-exp(-alpha_0-mu_0))*alpha_0/(alpha_0+mu_0);
PAttc=1-exp(-alpha_0-mu_0);

L[1]=1;
InitialPreprandialkillingEfficacy[1]=0;
Pf[1]=Pfc;
PAtt[1]=PAttc;

for (j in 1:tr) {
  L[j+1]=(log(2)^(1/kappa[j+1]))/beta[j+1];
  Pf[j+1]=(1-exp(-alpha_i[j+1] -mu_i[j+1])) * alpha_i[j+1]/(alpha_i[j+1]+mu_i[j+1]);
  PAtt[j+1]=1-exp(-alpha_i[j+1] -mu_i[j+1]);
  InitialPreprandialkillingEfficacy[j+1]=1-(Pf[j+1]/Pfc);
  InitialRepellentEfficacy[j+1]=1-(PAtt[j+1]/PAttc);
}
}
