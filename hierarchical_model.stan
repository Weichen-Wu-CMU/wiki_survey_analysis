//
// This stan program is for estimating the 
// hierarchical Thurstone (Salganik-Levy) 
// model to analyze wiki-survey data
data {
  int<lower=0> V; // number of votes
  int<lower=0> K; // number of ideas
  int<lower=0> N; // size of theta_v
  int Ns[K+1]; //Ns[k] means the number of theta_jk in theta_v
  int left_indices[V];
  int right_indices[V];
  int<lower = 0, upper = 1> y[V];
  real mu0[K];
  real tau0[K];
  real sigma;
}

// The parameters accepted by the model. 
parameters {
  real theta_v[N];
  real mu[K];
}

// The model to be estimated. 
model {
  for(k in 1:K)
  {
    mu[k] ~ normal(mu0[k],tau0[k]);
    for(i in Ns[k]:(Ns[k+1]-1))
    {
      theta_v[i] ~ normal(mu[k],sigma);
    }
  }
  for(i in 1:V)
  {
    y[i] ~ bernoulli(Phi_approx(theta_v[left_indices[i]] - theta_v[right_indices[i]]));
  }
}

//generated quantities: used for model evaluation
generated quantities {
  real deviance;
  real deviance_rep;
  int y_rep = 0;
  deviance = 0;
  deviance_rep = 0;
  for(i in 1:V)
  {
    deviance = deviance - 2 * bernoulli_lpmf(y[i]|
    Phi_approx(theta_v[left_indices[i]] - theta_v[right_indices[i]]));
    y_rep = bernoulli_rng(Phi_approx(theta_v[left_indices[i]] - theta_v[right_indices[i]]));
    deviance_rep = deviance_rep - 2 * bernoulli_lpmf(y_rep|
    Phi_approx(theta_v[left_indices[i]] - theta_v[right_indices[i]]));
  }
}
