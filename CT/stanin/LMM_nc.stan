// LMM with individuals means for each component
// as well as REs for subject and bigrams and by-subjects slopes
// Next model needs variance for each component

data {
	int<lower=1> N;                    // Number of observations
	real y[N];  		            //outcome
	int<lower=1> K; // N of components
	int<lower=1, upper=K> components[N];  //predictor

	int<lower=1> S;                  //number of subjects
	int<lower=1, upper=S> subj[N];   //subject id
	int<lower=1> B;                  //number of bigrams
	int<lower=1, upper=B> bigram[N];   //bigramid
}

parameters {
	real<lower=0> sigma;		// residual sd
	
  // For random effects
	vector[S] u; //subj intercepts
  real<lower=0> sigma_u;//subj sd
  
  vector[B] w; //bigram intercepts
  real<lower=0> sigma_w;//bigram sd

	// Parameters for non-centering
	real<lower = 0> mu_beta;
  real<lower = 0> sigma_beta;	
	vector[K] beta_raw;			// distributions
  
}

transformed parameters{
  vector[N] mu;
	vector[K] beta = mu_beta + sigma_beta * beta_raw;

  for(n in 1:N){
    mu[n] = beta[components[n]] + u[subj[n]] + w[bigram[n]];
  }
}

model {
  // Priors
  mu_beta ~ normal(5, 2.5);
  sigma_beta ~ normal(0, 2.5);
  beta_raw ~ normal(0, 1);
	sigma ~ normal(0, 1.5);
	
	// REs priors
  sigma_u ~ normal(0,2.5);
  u ~ normal(0, sigma_u); //subj random effects

  sigma_w ~ normal(0,2.5);
  w ~ normal(0, sigma_w); //bigram random effects

  // likelihood
  y ~ lognormal(mu, sigma);
}

generated quantities{
	real log_lik[N];
	vector[N] y_tilde;
	for (n in 1:N){
		log_lik[n] = lognormal_lpdf(y[n] | mu[n], sigma);
		y_tilde[n] = lognormal_rng(mu[n], sigma);
	}
}
