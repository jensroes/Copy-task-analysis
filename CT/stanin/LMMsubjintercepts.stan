// Null model with random intercepts for subjects and items
// weakly regulating priors

data {
	int<lower=1> N;                    // Number of observations
	real y[N];  		            //outcome
	
	int<lower=1> S;                  //number of subjects
	int<lower=1, upper=S> subj[N];   //subject id
	int<lower=1> B;                  //number of bigrams
	int<lower=1, upper=B> bigram[N];   //bigramid
}

parameters {
	real<lower=0> beta;			// intercept and slopes
	real<lower=0> sigma;		// residual sd
	
	  // For random effects
  vector[S] u; //subject intercepts
  real<lower=0> sigma_u;//subj sd
  
  vector[B] w; //bigram intercepts
  real<lower=0> sigma_w;//bigram sd

  
}

transformed parameters{
  vector[N] mu;
  for(n in 1:N){
    mu[n] = beta + u[subj[n]] + w[bigram[n]];
  }
}

model {
  // Priors
  beta ~ normal(5, 1);
	sigma ~ normal(0, 1.5);
	
	// REs priors
  sigma_u ~ normal(0, 2.5);
  u ~ normal(0, sigma_u); //subj random effects

  sigma_w ~ normal(0, 2.5);
  w ~ normal(0, sigma_w); //bigram random effects

  // likelihood
  y ~ lognormal(mu, sigma);
}

generated quantities{
	real log_lik[N];
	vector[N] y_tilde;
	vector[N] y_hat = mu;
	for (n in 1:N){
		log_lik[n] = lognormal_lpdf(y[n] | mu[n], sigma);
		y_tilde[n] = lognormal_rng(mu[n], sigma);
	}
}
