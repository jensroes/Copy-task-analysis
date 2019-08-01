/*
 Mixture model with 2 mixture components
 This model is adapted from Vasishth et al. 2007 
 Random intercepts and by-subject slopes adjustments for copy task component
 Random intercepts for bigrams
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
	real<lower=0> sigma[K];		// residual sd

	real<lower=0> sigma_diff[K];
	simplex[2] theta[K]; //probability of extreme values
  //real<lower=0.001, upper=0.999> theta[K];
	
  // For random effects
	vector<lower=0>[K] sigma_u;	// subj sd
	cholesky_factor_corr[K] L_u;
	matrix[K,S] z_u;
	
	vector[B] w; //bigrams intercepts
  real<lower=0> sigma_w;//bigram sd
	
}


transformed parameters{
  real mu1[N];
  real mu2[N];
  real sigmap_e_comp[N];
  real sigma_e_comp[N];
	real sigmap_e[K];
	real sigma_e[K];
	matrix[K,2] log_theta;
	matrix[N,2] log_theta_comp;
	matrix[S,K] u = (diag_pre_multiply(sigma_u, L_u) * z_u)';	// subj random effects

  for(k in 1:K){
  	log_theta[k,1] = log(theta[k,1]);
	  log_theta[k,2] = log1m(theta[k,1]);
	  sigmap_e[k] = sigma[k] + sigma_diff[k];
    sigma_e[k] = sigma[k] - sigma_diff[k];
  }
  
  for(n in 1:N){
    mu1[n] = beta[components[n]] + u[subj[n], components[n]] + w[bigram[n]];
    mu2[n] = beta[components[n]] + delta[components[n]] + u[subj[n], components[n]] + w[bigram[n]];
    for(i in 1:2){
      log_theta_comp[n,i] = log_theta[components[n],i]; 
    }
    sigmap_e_comp[n] = sigmap_e[components[n]];
    sigma_e_comp[n] = sigma_e[components[n]];
  }
}

model {
  vector[2] lp_parts[N];

  // Priors
  beta ~ normal(5, 1);
  delta ~ normal(0, 5);
  sigma ~ normal(0, 1.5);
  sigma_diff ~ normal(0, 1);

  for(k in 1:K){
    theta[k] ~ beta(2,2);
  }
  
	// REs priors
	sigma_u ~ normal(0,2.5);
	L_u ~ lkj_corr_cholesky(2);
	to_vector(z_u) ~ normal(0,1);	
	
	sigma_w ~ normal(0,2.5);
  w ~ normal(0, sigma_w); //bigram random effects
  
  // Likelihood	
  for (n in 1:N){
    lp_parts[n, 1] = log_theta_comp[n,1] + 
	      lognormal_lpdf(y[n] | mu1[n], sigma_e_comp[n]);
    lp_parts[n, 2] = log_theta_comp[n,2] + 
	      lognormal_lpdf(y[n] | mu2[n], sigmap_e_comp[n]);
    target += log_sum_exp(lp_parts[n]); 
  }
}

generated quantities{
	real log_lik[N];
	vector[N] y_tilde;
  real<lower=0,upper=1> theta_tilde; 
  real beta2[K];
  
  for(k in 1:K){
    beta2[k] = beta[k] + delta[k];
  }
  
  // likelihood: 
  for(n in 1:N){ 
	    log_lik[n] = log_sum_exp(
 		      log_theta_comp[n,1] + 
	    lognormal_lpdf(y[n] | mu1[n], sigma_e_comp[n]), 
 		      log_theta_comp[n,2] + 
	    lognormal_lpdf(y[n] | mu2[n], sigmap_e_comp[n])
	    );
       theta_tilde = bernoulli_rng(theta[components[n],1]); 
          if(theta_tilde) { 
              y_tilde[n] = lognormal_rng(mu1[n], sigma_e_comp[n]);
          }
          else{
              y_tilde[n] = lognormal_rng(mu2[n], sigmap_e_comp[n]);
          }
    }
}


