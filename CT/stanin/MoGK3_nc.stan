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
	vector<lower=0>[K] delta;			// distribution + extra component
	vector<lower=0>[K] gamma;			// distribution - extra component
	vector<lower=0>[K] sigma;		// residual sd
	real<lower=0> sigma_bar;
	real<lower=0> sigma_sigma;

	vector<lower=0>[K] sigma_diff;
	vector<lower=0>[K] sigma_diff2;

//	real<lower=0.001, upper=0.999> theta[3,K]; //probability of extreme values
//  matrix<lower=0.001, upper=0.999>[3, K] theta; //probability of extreme values
	simplex[3] theta[K]; //probability of extreme values
			
  // For random effects
  vector[S] u;
	real<lower=0> sigma_u;	// subj sd

	vector[B] w; //bigrams intercepts
  real<lower=0> sigma_w;//bigram sd
  
	// Parameters for non-centering
	real<lower = 0> mu_beta;
  real<lower = 0> sigma_beta;	
	vector[K] beta_raw;			// distributions
}


transformed parameters{
	real mu1[N];
	real mu2[N];
	real mu3[N];

	vector[K] sigmap_ee = sigma - sigma_diff2;
	vector[K] sigmap_e = sigma + sigma_diff;
	vector[K] sigma_e = sigma - sigma_diff;

  real sigmap_e_comp[N];
  real sigmap_ee_comp[N];
  real sigma_e_comp[N];
  
	matrix[N,3] log_theta_comp;
	matrix[K,3] log_theta;
	
	vector[K] beta = mu_beta + sigma_beta * beta_raw;


	for(k in 1:K){
	  for(i in 1:3){
  	  log_theta[k,i] =  log(theta[k,i]);
	  }
  }
  
  for(n in 1:N){
    mu1[n] = beta[components[n]] - gamma[components[n]] + u[subj[n]] + w[bigram[n]];
    mu2[n] = beta[components[n]] + u[subj[n]] + w[bigram[n]];
    mu3[n] = beta[components[n]] + delta[components[n]] + u[subj[n]] + w[bigram[n]];
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
  mu_beta ~ normal(5, 2.5);
  sigma_beta ~ normal(0, 2.5);
  beta_raw ~ normal(0, 1);

  delta ~ normal(0, 5);
  sigma_diff ~ normal(0, 1);
  sigma_diff2 ~ normal(0, 1);

  sigma_bar ~ normal(0, 1.5);
  sigma_sigma ~ normal(0, 1);
	sigma ~ normal(sigma_bar, sigma_sigma);

  for(k in 1:K){
    gamma[k] ~ uniform(0, beta[k]);
    theta[k] ~ dirichlet(rep_vector(2, 3));//beta(2, 2);
  }


	// REs priors
	sigma_u ~ normal(0,2.5);
  u ~ normal(0, sigma_u); //subj random effects

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
  vector[K] beta3 = beta + delta;
  vector[K] beta1 = beta - gamma;

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
