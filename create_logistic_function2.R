make_logistic2  <- function(max_VE,
                           time_to_max,
                           vacc_delay_mean,
                           vacc_delay_sd,
                           cluster_period = 30,
                           n_draws = 1) {
  
  # Delay in vaccination 
  
  #Draw this delay from a gamma distribution
  shape_vacc_delay = vacc_delay_mean^2/vacc_delay_sd^2
  scale_vacc_delay = vacc_delay_sd^2/vacc_delay_mean
  
  vacc_delay <- rgamma(n=n_draws, shape=shape_vacc_delay,  scale = scale_vacc_delay)
  
  day_delay_det <- round(vacc_delay)
  days_delay <- seq(from=1,to=day_delay_det)
  
  #Vaccine efficacy on these first few days is 0
  vacc_eff_start <- rep(0.000,day_delay_det)
  
  
  #Delay for vaccine to reach full efficacy
  delay_vacc_eff <- time_to_max
  days_delay_vacc_eff <- seq(from = 1, to = delay_vacc_eff) +
    day_delay_det
  
  
  
  #Maximum vaccine efficacy
  vacc_eff_max <- max_VE
  
  
  
  #Vaccine efficacy between 0 and maximum efficacy
  vac_eff_2 <- seq(from = 0, to = vacc_eff_max,length.out =delay_vacc_eff)
  
  #days at maximum 
  days_at_max <- max(0,
                     cluster_period - (day_delay_det + delay_vacc_eff) )
  
  
  if(days_at_max > 0){
  days_max_eff <- seq(from = 1, to = days_at_max) +
    day_delay_det + 
    delay_vacc_eff
  }  else{days_max_eff <- vector() }
  
  #vaccine efficacy during this period
  vacc_eff_max_period <- rep(vacc_eff_max,days_at_max)
  
  # Combine into dataframe
  
  #First the days
  day <- c(days_delay,  days_delay_vacc_eff, days_max_eff)
  
  if(length(day)!= cluster_period){
    day <- day[1:cluster_period]
    
  }
  
  #day <- c(days_delay, days_max_eff)
  
  #Then the vaccine efficacy on each day
  vaccine_efficacy <- c(vacc_eff_start,vac_eff_2,vacc_eff_max_period)
  
  if( length(vaccine_efficacy) != cluster_period){
    vaccine_efficacy <-  vaccine_efficacy[1:cluster_period]
    
  }
  
  vacc_eff_df <- data.frame(day,vaccine_efficacy)
  
  
  source("logistic_function2.R")
  
  log_params <- log.fit2('vaccine_efficacy','day',vacc_eff_df)
  
  phi1 <- log_params$`Logistic Curve`[1]
  phi2 <- log_params$`Logistic Curve`[2]
  phi3 <- log_params$`Logistic Curve`[3]
  
  #Find the average value
  
  VE_values <- vector()
  
  for (i in seq(vacc_eff_df$day)){
    VE_values[i] <- phi1/ (1 + exp(-(i - phi2)/phi3 ))
  }
  
  
  
  # Either obtain the mean or draw values from the distribution
  #Ask John/Thibaut
  
  vacc_eff <- mean(VE_values)
  
  #vacc_eff <-sample(VE_values, n_draws) 
  
  
  return(vacc_eff)
  
  
}


