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
	real<lower=0> beta[K];			// intercept
	real<lower=0> sigma[K];		// residual sd
	
  // For random effects
	vector<lower=0>[K] sigma_u;	// subj sd
	cholesky_factor_corr[K] L_u;
	matrix[K,S] z_u; // tmp for subject intercepts and slopes
  
  vector[B] w; //bigram intercepts
  real<lower=0> sigma_w;//bigram sd

  
}

transformed parameters{
  vector[N] mu;
  vector[N] sigma_comp;
  matrix[S,K] u = (diag_pre_multiply(sigma_u, L_u) * z_u)';	// subj random effects

  for(n in 1:N){
    mu[n] = beta[components[n]] + u[subj[n], components[n]] + w[bigram[n]];
    sigma_comp[n] = sigma[components[n]];
  }
}

model {
  // Priors
  beta ~ normal(5, 1);
	sigma ~ normal(0, 1.5);
	
	// REs priors
	sigma_u ~ normal(0,2.5);
	L_u ~ lkj_corr_cholesky(3);
	to_vector(z_u) ~ normal(0,1);	

  sigma_w ~ normal(0,2.5);
  w ~ normal(0, sigma_w); //bigram random effects

  // likelihood
   y ~ lognormal(mu, sigma_comp);
}

generated quantities{
	real log_lik[N];
	vector[N] y_tilde;
	for (n in 1:N){
		log_lik[n] = lognormal_lpdf(y[n] | beta[components[n]] + u[subj[n], components[n]] + w[bigram[n]], sigma[components[n]]);
		y_tilde[n] = lognormal_rng(beta[components[n]] + u[subj[n], components[n]] + w[bigram[n]], sigma[components[n]]);
	}
}
