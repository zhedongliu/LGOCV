data {
  int<lower=0> n; 
  int<lower=0> m;
  real y[n];
  int id[n];
}
parameters {
  real mu_group[m];
  real mu;
  real logprec_group;
}

transformed parameters {
  real sigma_group; 
  real prec;
  prec = exp(logprec_group);
  sigma_group = sqrt(1/prec);
}

model {
  logprec_group ~ normal(0,100);
  mu_group ~ normal(0,sigma_group);
  mu ~ normal(0,100);
  for(i in 1:n){
    y[i] ~ normal(mu+mu_group[id[i]],.1);
  }
}

