/*
 Mixture model with 3 mixture components 
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
	real<lower=0> gamma[K];			// distribution - extra component
	real<lower=0> sigma[K];		// residual sd
	real<lower=0> sigma_diff[K];
	real<lower=0> sigma_diff2[K];

//	real<lower=0.001, upper=0.999> theta[3,K]; //probability of extreme values
//  matrix<lower=0.001, upper=0.999>[3, K] theta; //probability of extreme values
	simplex[3] theta[K]; //probability of extreme values
			
  // For random effects
	vector<lower=0>[K] sigma_u;	// subj sd
	cholesky_factor_corr[K] L_u;
	matrix[K,S] z_u;

	vector[B] w; //bigrams intercepts
  real<lower=0> sigma_w;//bigram sd
}


transformed parameters{
	real sigmap_ee[K];
	real sigmap_e[K];
	real sigma_e[K];
	real mu1[N];
	real mu2[N];
	real mu3[N];
  real sigmap_e_comp[N];
  real sigmap_ee_comp[N];
  real sigma_e_comp[N];
	matrix[N,3] log_theta_comp;
	matrix[K,3] log_theta;
	matrix[S,K] u = (diag_pre_multiply(sigma_u, L_u) * z_u)';	// subj random effects


	for(k in 1:K){
	  for(i in 1:3){
  	  log_theta[k,i] =  log(theta[k,i]);
	  }
    sigmap_e[k] = sigma[k] + sigma_diff[k];
    sigma_e[k] = sigma[k] - sigma_diff[k];
    sigmap_ee[k] = sigma[k] - sigma_diff2[k];
  }
  
  for(n in 1:N){
    mu1[n] = beta[components[n]] - gamma[components[n]] + u[subj[n], components[n]] + w[bigram[n]];
    mu2[n] = beta[components[n]] + u[subj[n], components[n]] + w[bigram[n]];
    mu3[n] = beta[components[n]] + delta[components[n]] + u[subj[n], components[n]] + w[bigram[n]];
    sigmap_ee_comp[n] = sigmap_ee[components[n]];
    sigmap_e_comp[n] = sigmap_e[components[n]];
    sigma_e_comp[n] = sigma_e[components[n]];
    for(i in 1:3){
      log_theta_comp[n,i] = log_theta[components[n],i];
    }
  }
}

model {
  vector[3] lp_parts[N];

  // Priors
	beta ~ normal(5, 1);
  delta ~ normal(0, 5);
  sigma_diff ~ normal(0, 1);
  sigma_diff2 ~ normal(0, 1);
  sigma ~ normal(0, 1.5);

  for(k in 1:K){
    gamma[k] ~ uniform(0, beta[k]);
    theta[k] ~ dirichlet(rep_vector(2, 3));//beta(2, 2);
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
  	    lognormal_lpdf(y[n] | mu1[n], sigmap_ee_comp[n]);
  	  lp_parts[n, 2] = log_theta_comp[n, 2] + 
  	    lognormal_lpdf(y[n] | mu2[n], sigma_e_comp[n]);
  	  lp_parts[n, 3] = log_theta_comp[n, 3] + 
  	    lognormal_lpdf(y[n] | mu3[n], sigmap_e_comp[n]);
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
      ps[1] = log_theta_comp[n,1] + lognormal_lpdf(y[n] | mu1[n], sigmap_ee_comp[n]);
      ps[2] = log_theta_comp[n,2] + lognormal_lpdf(y[n] | mu2[n], sigma_e_comp[n]);
      ps[3] = log_theta_comp[n,3] + lognormal_lpdf(y[n] | mu3[n], sigmap_e_comp[n]);
      log_lik[n] = log_sum_exp(ps);
  	  theta_tilde = categorical_rng(theta[components[n]]); 
      if(theta_tilde == 1){
        y_tilde[n] = lognormal_rng(mu1[n], sigmap_ee_comp[n]);
      }
      if(theta_tilde == 2){
        y_tilde[n] = lognormal_rng(mu2[n], sigma_e_comp[n]);
      }
      if(theta_tilde == 3){
        y_tilde[n] = lognormal_rng(mu3[n], sigmap_e_comp[n]);
      }
  	}
}
