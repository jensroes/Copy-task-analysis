/*
 Mixture model with 3 mixture components and random intercepts and slopes (for subject)
 This model is adapted from Vasishth et al. 2007 
*/

data {
	int<lower=1> N;                    // Number of observations
	real y[N];  		            //outcome
	int<lower=1> K;             // N of copy task components 
	int<lower=1, upper=K> components[N];  //predictor
	
	int<lower=1> S;                  //number of subjects
	int<lower=1, upper=S> subj[N];   //subject id

	int<lower=1> B;                  //number of bigrams
	int<lower=1, upper=B> bigram[N];   //bigram id
}


parameters {
	real<lower=0> beta[K];			// distributions
	real<lower=0> delta[K];			// distribution + extra component
	real<lower=0> gamma[K];			// distribution - extra component
	real<lower=0> sigma[K,3];		// residual sd

	simplex[3] theta[K]; //probability of extreme values
			
  // For random effects
	vector<lower=0>[K] sigma_u;	// subj sd
	cholesky_factor_corr[K] L_u;
	matrix[K,S] z_u;
	vector[B] w; //bigrams intercepts
  real<lower=0> sigma_w;//bigram sd
}


transformed parameters{
	real<lower=0> mu1[N];
	real<lower=0> mu2[N];
	real<lower=0> mu3[N];

  matrix[N,3] sigma_comp;
	matrix[N,3] log_theta_comp;
	matrix[K,3] log_theta;
	matrix[S,K] u = (diag_pre_multiply(sigma_u, L_u) * z_u)';	// subj random effects


	for(k in 1:K){
	  for(i in 1:3){
  	  log_theta[k,i] =  log(theta[k,i]);
	  }
  }
  
  for(n in 1:N){
    mu1[n] = beta[components[n]] - gamma[components[n]] + u[subj[n], components[n]] + w[bigram[n]];
    mu2[n] = beta[components[n]] + u[subj[n], components[n]] + w[bigram[n]];
    mu3[n] = beta[components[n]] + delta[components[n]] + u[subj[n], components[n]] + w[bigram[n]];
    for(i in 1:3){
      sigma_comp[n,i] = sigma[components[n],i];
      log_theta_comp[n,i] = log_theta[components[n],i];
    }
  }
}

model {
  vector[3] lp_parts[N];

  // Priors
  for(k in 1:K){
  	beta[k] ~ normal(5, 5);
    delta[k] ~ normal(0, 10);
    gamma[k] ~ uniform(0, beta[k]);
    for(i in 1:3){
      sigma[k,i] ~ normal(0, 10);
    }
    theta[k] ~ dirichlet(rep_vector(2,3));//beta(2, 2);
  }

	// REs priors
	sigma_u ~ normal(0,2.5);
	L_u ~ lkj_corr_cholesky(2);
	to_vector(z_u) ~ normal(0,1);	 //subj random effects
	
	sigma_w ~ normal(0,2.5);
  w ~ normal(0, sigma_w); //bigram random effects

  // Likelihood	
	for(n in 1:N){
  	  lp_parts[n, 1] = log_theta_comp[n, 1] + 
  	    lognormal_lpdf(y[n] | mu1[n], sigma_comp[n,1]);
  	  lp_parts[n, 2] = log_theta_comp[n, 2] + 
  	    lognormal_lpdf(y[n] | mu2[n], sigma_comp[n,2]);
  	  lp_parts[n, 3] = log_theta_comp[n, 3] + 
  	    lognormal_lpdf(y[n] | mu3[n], sigma_comp[n,3]);
  	  target += log_sum_exp(lp_parts[n]);
  }
}

generated quantities{
	real log_lik[N];
	vector[N] y_tilde;
	real ps[3];
	
  real<lower=1,upper=3> theta_tilde; 
  real beta3[K];
  real beta1[K];
  
  for(k in 1:K){
    beta3[k] = beta[k] + delta[k];
    beta1[k] = beta[k] - gamma[k];
  }
  
  // likelihood: 
  for(n in 1:N){ 
      ps[1] = log_theta_comp[n,1] + lognormal_lpdf(y[n] | mu1[n], sigma_comp[n,1]);
      ps[2] = log_theta_comp[n,2] + lognormal_lpdf(y[n] | mu2[n], sigma_comp[n,2]);
      ps[3] = log_theta_comp[n,3] + lognormal_lpdf(y[n] | mu3[n], sigma_comp[n,3]);
      log_lik[n] = log_sum_exp(ps);
  	  theta_tilde = categorical_rng(theta[components[n]]); 
      if(theta_tilde == 1){
        y_tilde[n] = lognormal_rng(mu1[n], sigma_comp[n,1]);
      }
      if(theta_tilde == 2){
        y_tilde[n] = lognormal_rng(mu2[n], sigma_comp[n,2]);
      }
      if(theta_tilde == 3){
        y_tilde[n] = lognormal_rng(mu3[n], sigma_comp[n,3]);
      }
  	}
}



