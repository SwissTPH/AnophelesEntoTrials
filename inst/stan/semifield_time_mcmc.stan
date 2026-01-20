// https://pmc.ncbi.nlm.nih.gov/articles/PMC7819244/#_ad93_
// Stan model for EACoMoPP semifield experiments, control 

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
int<lower=0> y0[n,K]; // observed outcome vectors (H1, ..., H4, T, L) for n experiments from control arm
int<lower=0> y1[n,K]; // observed outcome vectors (H1, ..., H4, T, L) for n experiments from intervention arm 
real<lower=0> priorsigma_mean_logrates; // sigma of the normal priors on a and m
real<lower=0> priorsigma_pikappa; // sigma of the priors on pi and kappa
real<lower=0> hierarchy; // variance for phi and psi, limit --> 0 means eliminating hierarchy
}
parameters{
  real<lower=0> rho; 
  real a;
  real b;
  real m;
  real<lower=0> sigma_a;
  real<lower=0> sigma_b;
  real<lower=0> sigma_m;
  vector[n] phi;
  vector[n] eta;
  vector[n] psi;
}
transformed parameters{
  // declare: control rates
  vector<lower=0>[n] alpha_H_k0;
  vector<lower=0>[n] alpha_T_k0;
  vector<lower=0>[n] mu_k0;
  // declare: intervention rates
  vector<lower=0>[n] alpha_H_k1;
  vector<lower=0>[n] alpha_T_k1;
  vector<lower=0>[n] mu_k1;
  // declare: probabilities of multinomial model for semifield experiment
  matrix[6,n] theta0; // be carefull: one column is one replicate!!!
  matrix[6,n] theta1; // be carefull: one column is one replicate!!!
  // define: control rates
  alpha_H_k0 = exp(a + phi *sigma_a);
  alpha_T_k0 = exp(b + eta *sigma_b);
  mu_k0 = exp(m + psi *sigma_m);
  // define: intervention rates
  alpha_H_k1 = alpha_H_k0;
  alpha_T_k1 = rho * alpha_H_k0;
  mu_k1 = mu_k0;
  // define probabilities of multinomial model for semifield experiment of control arm
    // be careful: one column is one replicate!!!
  theta0[1,] = (p_H1(alpha_H_k0, alpha_T_k0, mu_k0))';
  theta0[2,] = (p_H2(alpha_H_k0, alpha_T_k0, mu_k0))';
  theta0[3,] = (p_H3(alpha_H_k0, alpha_T_k0, mu_k0))';
  theta0[4,] = (p_H4(alpha_H_k0, alpha_T_k0, mu_k0))';
  theta0[5,] = (p_T(alpha_H_k0, alpha_T_k0, mu_k0))';
  theta0[6,] = (p_L(alpha_H_k0, alpha_T_k0, mu_k0))';
  // define probabilities of multinomial model for semifield experiment of intervention arm
    // be careful: one column is one replicate!!!
  theta1[1,] = (p_H1(alpha_H_k1, alpha_T_k1, mu_k1))';
  theta1[2,] = (p_H2(alpha_H_k1, alpha_T_k1, mu_k1))';
  theta1[3,] = (p_H3(alpha_H_k1, alpha_T_k1, mu_k1))';
  theta1[4,] = (p_H4(alpha_H_k1, alpha_T_k1, mu_k1))';
  theta1[5,] = (p_T(alpha_H_k1, alpha_T_k1, mu_k1))';
  theta1[6,] = (p_L(alpha_H_k1, alpha_T_k1, mu_k1))';
}
model{
//declare local variables
  vector[6] theta0_i;
  vector[6] theta1_i;
//priors
  // priors on effect parameters
  target += lognormal_lpdf(rho | 0, priorsigma_pikappa);
  // priors on means of logrates
  a ~ normal(0,priorsigma_mean_logrates);  
  b ~ normal(0,priorsigma_mean_logrates); 
  m ~ normal(0,priorsigma_mean_logrates);  
  // priors on phi and psi (normalised deviation from mean of logrates for each group)
  phi ~ normal(0,1);
  eta ~ normal(0,1);
  psi ~ normal(0,1);
  // priors on sigma_a and sigma_m (variation of deviation from mean of logrates for each group)
  target += 2*cauchy_lpdf(sigma_a | 0,hierarchy);
  target += 2*cauchy_lpdf(sigma_b | 0,hierarchy);
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
  real<lower=0> alpha_T0; // 
  real<lower=0> mu0; //
  real<lower=0> alpha_H1; // 
  real<lower=0> alpha_T1; // 
  real<lower=0> mu1; //
    //transform: means of rates
  alpha_H0 = exp(a + sigma_a^2/2);
  alpha_T0 = exp(b + sigma_b^2/2);
  mu0 = exp(m + sigma_m^2/2);
  alpha_H1 = alpha_H0;
  alpha_T1 = rho * alpha_H0;
  mu1 = mu0;
}
