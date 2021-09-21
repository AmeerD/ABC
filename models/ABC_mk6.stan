functions {
  vector rescale_bias(vector x){
    int len = num_elements(x);
    real c = inv(len);
    vector[len + 1] beta_full;
    matrix[len + 1, len] trans= rep_matrix(-c, len + 1, len);
    
    for (i in 1:len){
      trans[i,i] = 1 - (i-1)*c;
      trans[i, i+1:len] = rep_row_vector(0, len - i);
    }
    
    beta_full = trans * x;
    
    // beta_full[1] = x[1];
    // for (i in 2:len){
    //   row_vector[len] temp = rep_row_vector(0, len);
    //   temp[1:(i-1)] = rep_row_vector(-c, i-1);
    //   temp[i] = 1 - (i-1)*c;
    //   beta_full[i] = temp * x;
    // }
    // beta_full[len + 1] = rep_row_vector(-c, len) * x;
    
    return beta_full;
  }
}

data {
  int n;
  int nyears;
  int n_country;
  int n_survey;
  int survey_count[n_country]; //Vector of country identifiers for each survey
  int year[n];
  int country[n];  
  int survey[n];
  vector[n] obsage;
  int recondist[n];
  int truage5mlt[n];
  vector[n] value;
  vector[n] se_q2;
  int haslate2[n_country];
  int hasvlate2[n_country];
  vector[n_country] average;
}

transformed data {
  int agg_survey_count[n_country + 1];
  int recondist3[n];
  vector[n] reconvar;
  vector[n] halfobsage;
  int nlate = sum(haslate2);
  int nvlate = sum(hasvlate2);
  int lateidx[n_country];
  int vlateidx[n_country];
  int counter = 1;
  int vcounter = 1;
  
  agg_survey_count[1] = 0;
  for (i in 1:n_country){
    agg_survey_count[i+1] = sum(survey_count[1:i]);
  }
  
  for (i in 1:n){
    recondist3[i] = min(3, recondist[i]);
    reconvar[i] = 1 + 0.05*recondist[i];
    halfobsage[i] = obsage[i]/2;
  }
  
  for (i in 1:n_country){
    if (haslate2[i] == 0){
      lateidx[i] = nlate + 1;
    } else {
      lateidx[i] = counter;
      counter += 1;
    }
    if (hasvlate2[i] == 0){
      vlateidx[i] = nvlate + 1;
    } else {
      vlateidx[i] = vcounter;
      vcounter += 1;
    }
  }
  
}

parameters {
  matrix[n_country, nyears] e_ct; //AR(1) white noise
  vector[n_country] muzero; //intercept
  vector<lower=0>[n_country] drift;  //latent drift
  //vector<lower=0, upper=1>[n_country] rho_beta; //vector<lower=0, upper=1>[n_country] rho_beta; //autocorrelation beta term
  
  vector<lower=0, upper=1>[n_country] rho;
  
  vector<lower=-pi()/2, upper=pi()/2>[n_survey - n_country] beta_s_unif;
  vector<lower=0>[nlate] late_short;
  vector<lower=0>[nvlate] vlate_short;
  vector<lower=0>[n_country] mult5err; //Lower bound to correct for negative effects of the 5-year error
  vector<lower=0>[n_survey] sigma_s;
  //real mu_ln;
  //real<lower=0> sigma_ln;
  real<lower=0> sigma_beta_s;
  real<lower=0> sigma_e;
  real<lower=0> sigma_late;
  real<lower=0> sigma_vlate;
  real<lower=0> lambda_5;

  real alpha_d;
  real<lower=0> beta_d;
  real<lower=0> tau;
  vector<lower=0>[n_country] lambda_rho;
}

transformed parameters {
  matrix[n_country, nyears] mu_ct; 
  vector[n_survey - n_country] beta_s_short;
  vector[n_survey] beta_s = rep_vector(0, n_survey);
  //vector<lower=-1, upper=1>[n_country] rho = 2*(rho_beta - 0.5); //autocorrelation term
  vector<lower=0>[nlate+1] late_short2;
  vector<lower=0>[nvlate+1] vlate_short2;
  
  mu_ct[1:n_country, 1] = muzero + e_ct[1:n_country, 1] + drift;
  mu_ct[1:n_country, 2] = mu_ct[1:n_country, 1] + e_ct[1:n_country, 2] + rho .* (mu_ct[1:n_country, 1] - muzero) + drift;
  for (i in 3:nyears){
    mu_ct[1:n_country, i] = mu_ct[1:n_country, i-1] + rho .* (mu_ct[1:n_country, i-1] - mu_ct[1:n_country, i-2]) + e_ct[1:n_country, i] + drift;
  }
  
  beta_s_short = tan(beta_s_unif) * sigma_beta_s; //cauchy(0,sigma_beta_s)

  for (i in 1:n_country){
    if (survey_count[i] != 1){
      beta_s[(agg_survey_count[i]+1):agg_survey_count[i+1]] = rescale_bias(beta_s_short[(agg_survey_count[i]-i+2):(agg_survey_count[i+1]-i)]);
    }
  }
  
  late_short2[1:nlate] = late_short;
  late_short2[nlate+1] = 0;
  vlate_short2[1:nvlate] = vlate_short;
  vlate_short2[nvlate+1] = 0;
}

model {

  //Prior Distributions
  sigma_e ~ gamma(2, 0.1);
  to_vector(e_ct) ~ normal(0, sigma_e); 
  muzero ~ normal(0, 10);
  
  tau ~ normal(0,0.01);
  lambda_rho ~ student_t(4, 0, 1);
  
  rho ~ normal(0, tau * lambda_rho);
  
  drift ~ lognormal(alpha_d, beta_d);
  alpha_d ~ normal(0, 1);
  beta_d ~ normal(0, 1);
  
  beta_s_unif ~ uniform(-pi()/2, pi()/2);
  sigma_beta_s ~ normal(0, 0.25);
  sigma_s ~ gamma(2, 4); 
  
  sigma_late ~ std_normal();
  sigma_vlate ~ std_normal();
  late_short ~ normal(0, sigma_late); 
  vlate_short ~ normal(0, sigma_vlate); 
  
  lambda_5 ~ normal(0, 50);
  mult5err ~ exponential(lambda_5);
  
  //Data Generation
  {
  vector[n] yhat;
  vector[n] sigma_total;
  
  for(i in 1:n){
    yhat[i] = mu_ct[country[i], year[i]] + beta_s[survey[i]] 
    + vlate_short2[vlateidx[country[i]]] * recondist3[i] - late_short2[lateidx[country[i]]] .* obsage[i]
    - mult5err[country[i]] * truage5mlt[i];
    
    sigma_total[i] = sqrt(se_q2[i]^2 + reconvar[i]*sigma_s[survey[i]]^2);
  }
  
  value ~ normal(yhat, sigma_total);
  }
  
}

generated quantities {
  vector[n_country] late;
  vector[n_country] vlate;
  vector[n] yhat;
  vector[n] yrepl;
  vector[n] ll;
  
  for (i in 1:n_country){
    late[i] = late_short2[lateidx[i]];
    vlate[i] = vlate_short2[vlateidx[i]];
  }
  
  for(i in 1:n){
    real sd_total = sqrt(se_q2[i]^2 + reconvar[i]*sigma_s[survey[i]]^2);
    
    yhat[i] = mu_ct[country[i], year[i]] + beta_s[survey[i]] 
    + vlate_short2[vlateidx[country[i]]] * recondist3[i] - late_short2[lateidx[country[i]]] .* obsage[i]
    - mult5err[country[i]] * truage5mlt[i];
    
    yrepl[i] = normal_rng(yhat[i], sd_total);
    ll[i] = normal_lpdf(value[i] | yhat[i], sd_total);
  }
}
