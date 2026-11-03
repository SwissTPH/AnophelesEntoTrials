// Denz et al. (2021) Parasites & Vectors 14:64, https://doi.org/10.1186/s13071-020-04560-x
// Stan model for EACoMoPP semifield experiments, matched control and intervention arms
// alpha_H in the article is alpha here and alpha_T in the article is beta here
// Trap availability (beta0, beta1) is set as data: 0 for a spatial repellent ("push") only
// InitialRepellencyRate is pi and KillingDuringHostSeeking is kappa in the article

functions{
  // probabilities in the continuous-time Markov chain model for mosquito foraging behaviour
  vector P_A(vector alpha, vector beta, vector mu, real t){
    return exp( -(alpha+beta+mu) *t );
  }
  vector P_H(vector alpha, vector beta, vector mu, real t){
    return (1-P_A(alpha, beta, mu, t)) .*alpha ./(alpha +beta +mu);
  }
  vector P_T(vector alpha, vector beta, vector mu, real t){
    return (1-P_A(alpha, beta, mu, t)) .*beta ./(alpha +beta +mu);
  }
  vector P_M(vector alpha, vector beta, vector mu, real t){
    return (1-P_A(alpha, beta, mu, t)) .*mu ./(alpha +beta +mu);
  }
  // probabilities of the multinomial model for semi-field experiments
  vector p_H1(vector alpha, vector beta, vector mu){
    return P_H(alpha, beta, mu, 1);
  }
  vector p_H2(vector alpha, vector beta, vector mu){
    return P_A(alpha, beta, mu, 1) .*P_H(alpha, beta, mu, 1);
  }
  vector p_H3(vector alpha, vector beta, vector mu){
    return P_A(alpha, beta, mu, 1) .*P_A(alpha, beta, mu, 1) .*P_H(alpha, beta, mu, 1);
  }
  vector p_H4(vector alpha, vector beta, vector mu){
    return P_A(alpha, beta, mu, 1) .*P_A(alpha, beta, mu, 1) .*P_A(alpha, beta, mu, 1) .*P_H(alpha, beta, mu, 9);
  }
  vector p_T(vector alpha, vector beta, vector mu){
    return P_T(alpha, beta, mu, 12);
  }
  vector p_L(vector alpha, vector beta, vector mu){
    return 1 -p_H1(alpha, beta, mu) -p_H2(alpha, beta, mu) -p_H3(alpha, beta, mu) -p_H4(alpha, beta, mu) -p_T(alpha, beta, mu);
  }
}
data{
int<lower=0> n; // number of experiments per arm
int<lower=0> K; // number of observed outcomes (H1, ..., H4, T, L) per experiment
array[n, K] int<lower=0> y0; // observed outcome vectors (H1, ..., H4, T, L) for n experiments from control arm
array[n, K] int<lower=0> y1; // observed outcome vectors (H1, ..., H4, T, L) for n experiments from intervention arm
vector<lower=0>[n] beta0; // trap availability set for controls arm
vector<lower=0>[n] beta1; // trap availability set for intervention arm
real<lower=0> priorsigma_mean_logrates; // sigma of the normal priors on a and m
real<lower=0> priorsigma_pikappa; // sigma of the priors on pi and kappa
real<lower=0> hierarchy; // variance for phi and psi, limit --> 0 means eliminating hierarchy
}
parameters{
  real<upper=1> InitialRepellencyRate;
  real<lower=0> KillingDuringHostSeeking; // forcing > 0 for now since otherwise problem!!
  real a;
  real m;
  real<lower=0> sigma_a;
  real<lower=0> sigma_m;
  vector[n] phi;
  vector[n] psi;
}
transformed parameters{
  // declare: control rates
  vector<lower=0>[n] alpha_H_k0;
  vector<lower=0>[n] mu_k0;
  // declare: intervention rates
  vector<lower=0>[n] alpha_H_k1;
  vector<lower=0>[n] mu_k1;
  // declare: probabilities of multinomial model for semifield experiment
  matrix[6,n] theta0; // be carefull: one column is one replicate!!!
  matrix[6,n] theta1; // be carefull: one column is one replicate!!!
  // transform: control rates
  alpha_H_k0 = exp(a + phi *sigma_a);
  mu_k0 = exp(m + psi *sigma_m);
  // transform and define: intervention rates
  alpha_H_k1 = (1 -InitialRepellencyRate) *alpha_H_k0;
  mu_k1 = mu_k0 + KillingDuringHostSeeking * alpha_H_k0;
  // define probabilities of multinomial model for semifield experiment of control arm
    // be careful: one column is one replicate!!!
  theta0[1,] = (p_H1(alpha_H_k0, beta0, mu_k0))';
  theta0[2,] = (p_H2(alpha_H_k0, beta0, mu_k0))';
  theta0[3,] = (p_H3(alpha_H_k0, beta0, mu_k0))';
  theta0[4,] = (p_H4(alpha_H_k0, beta0, mu_k0))';
  theta0[5,] = (p_T(alpha_H_k0, beta0, mu_k0))';
  theta0[6,] = (p_L(alpha_H_k0, beta0, mu_k0))';
  // define probabilities of multinomial model for semifield experiment of intervention arm
    // be careful: one column is one replicate!!!
  theta1[1,] = (p_H1(alpha_H_k1, beta1, mu_k1))';
  theta1[2,] = (p_H2(alpha_H_k1, beta1, mu_k1))';
  theta1[3,] = (p_H3(alpha_H_k1, beta1, mu_k1))';
  theta1[4,] = (p_H4(alpha_H_k1, beta1, mu_k1))';
  theta1[5,] = (p_T(alpha_H_k1, beta1, mu_k1))';
  theta1[6,] = (p_L(alpha_H_k1, beta1, mu_k1))';
}
model{
//declare local variables
  vector[6] theta0_i;
  vector[6] theta1_i;
//priors
  // priors on effect parameters
  target += lognormal_lpdf(1 -InitialRepellencyRate | 0, priorsigma_pikappa);
  target += 2* cauchy_lpdf(KillingDuringHostSeeking | 0,priorsigma_pikappa);
  // priors on means of logrates
  a ~ normal(0,priorsigma_mean_logrates);
  m ~ normal(0,priorsigma_mean_logrates);
  // priors on phi and psi (normalised deviation from mean of logrates for each group)
  phi ~ normal(0,1);
  psi ~ normal(0,1);
  // priors on sigma_a and sigma_m (variation of deviation from mean of logrates for each group)
  target += 2*cauchy_lpdf(sigma_a | 0,hierarchy);
  target += 2*cauchy_lpdf(sigma_m | 0,hierarchy);
 // loop for defining model
  for (i in 1:n){
    theta0_i = theta0[,i];
    theta1_i = theta1[,i];
    target += multinomial_lpmf( y0[i,] | theta0_i);
    target += multinomial_lpmf( y1[i,] | theta1_i);
  }
}
generated quantities {
  // declare: means of rates
  real<lower=0> alpha_H0; //
  real<lower=0>  mu0; //
  real<lower=0> alpha_H1; //
  real<lower=0>  mu1; //
    //transform: means of rates
  alpha_H0 = exp(a + sigma_a^2/2);
  mu0 = exp(m + sigma_m^2/2);
  alpha_H1 = (1 -InitialRepellencyRate) *alpha_H0;
  mu1 = mu0 + KillingDuringHostSeeking * alpha_H0;
}
