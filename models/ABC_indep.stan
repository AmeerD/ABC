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
  vector[n_country] drift;  //latent drift
  vector<lower=-pi()/2, upper=pi()/2>[n_survey - n_country] beta_s_unif;
  vector<lower=0>[nlate] late_short;
  vector<lower=0>[nvlate] vlate_short;
  vector<lower=0>[n_country] mult5err; //Lower bound to correct for negative effects of the 5-year error
  vector<lower=0>[n_survey] sigma_s;
  real<lower=0> sigma_beta_s;
  real<lower=0> sigma_e;
}

transformed parameters {
  matrix[n_country, nyears] mu_ct; 
  vector[n_survey - n_country] beta_s_short;
  vector[n_survey] beta_s;
  vector<lower=0>[nlate+1] late_short2;
  vector<lower=0>[nvlate+1] vlate_short2;
  
  mu_ct[1:n_country, 1] = muzero + e_ct[1:n_country, 1] + drift;
  for (i in 2:nyears){
    mu_ct[1:n_country, i] = mu_ct[1:n_country, i-1] + e_ct[1:n_country, i] + drift;
  }
  
  beta_s_short = tan(beta_s_unif) * sigma_beta_s; //cauchy(0,sigma_beta_s)

  for (i in 1:n_country){
    beta_s[(agg_survey_count[i]+1):agg_survey_count[i+1]-1] = beta_s_short[(agg_survey_count[i]-i+2):(agg_survey_count[i+1]-i)];
    beta_s[agg_survey_count[i+1]] = -sum(beta_s_short[(agg_survey_count[i]-i+2):(agg_survey_count[i+1]-i)]);
  }
  
  late_short2[1:nlate] = late_short;
  late_short2[nlate+1] = 0;
  vlate_short2[1:nvlate] = vlate_short;
  vlate_short2[nvlate+1] = 0;
}

model {

  //Prior Distributions
  sigma_e ~ gamma(2,0.1);
  to_vector(e_ct) ~ normal(0, sigma_e); 
  muzero ~ normal(0, 10);
  drift ~ exp_mod_normal(0, 0.025, 2);   
  beta_s_unif ~ uniform(-pi()/2, pi()/2);
  late_short ~ normal(0, 0.25); 
  vlate_short ~ normal(0, 0.25); 
  mult5err ~ exponential(10);
  sigma_beta_s ~ normal(0, 0.25);
  sigma_s ~ gamma(2, 4); 
  
  //Data Generation
  {
  vector[n] yhat;
  vector[n] sigma_total;
  
  
  for(i in 1:n){
    yhat[i] = mu_ct[country[i], year[i]] + beta_s[survey[i]] 
    + vlate_short2[vlateidx[country[i]]] * recondist3[i] - late_short2[lateidx[country[i]]] .* halfobsage[i]
    - mult5err[country[i]] * truage5mlt[i];
    
    sigma_total[i] = sqrt(se_q2[i]^2 + reconvar[i]*sigma_s[survey[i]]^2);
  }
  
  value ~ normal(yhat, sigma_total);
  }
  
}
generated quantities {
  vector[n_country] late;
  vector[n_country] vlate;
  //vector[n] yhat;
  //vector[n] yrepl;
  //vector[n] ll;
  
  for (i in 1:n_country){
    late[i] = late_short2[lateidx[i]];
    vlate[i] = vlate_short2[vlateidx[i]];
  }
  
  //for(i in 1:n){
  //  real sd_total = sqrt(se_q2[i]^2 + reconvar[i]*sigma_s[survey[i]]^2);
    
  //  yhat[i] = mu_ct[country[i], year[i]] + beta_s[survey[i]] 
  //  + vlate_short2[vlateidx[country[i]]] * recondist3[i] - late_short2[lateidx[country[i]]] .* halfobsage[i]
  //  - mult5err[country[i]] * truage5mlt[i];
    
  //  yrepl[i] = normal_rng(yhat[i], sd_total);
  //  ll[i] = normal_lpdf(value[i] | yhat[i], sd_total);
  //}
}
