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
	real<lower=0> beta;			// intercept
	real<lower=0> sigma;		// residual sd
	
  // For random effects
	vector<lower=0>[K] sigma_u;	// subj sd
	cholesky_factor_corr[K] L_u;
	matrix[K,S] z_u; // tmp for subject intercepts and slopes
  
  vector[B] w; //bigram intercepts
  real<lower=0> sigma_w;//bigram sd

  
}

transformed parameters{
  vector[N] mu;
	matrix[S,K] u;
	
	u = (diag_pre_multiply(sigma_u, L_u) * z_u)';	// subj random effects

  for(n in 1:N){
    mu[n] = beta + u[subj[n], components[n]] + w[bigram[n]];
  }
}

model {
  // Priors
  beta ~ normal(5, 1);
	sigma ~ cauchy(0, 2.5);
	
	// REs priors
	sigma_u ~ normal(0, 2.5);
	L_u ~ lkj_corr_cholesky(3);
	to_vector(z_u) ~ normal(0,1);	

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
