
# Toy example of a branching process simulator

# This script illustrates an iteration of a branching process simulator where
# individual features can be tracked for the cases, which impact their
# infectivity.

# Thibaut Jombart 2021 - MIT License


# ------------
# Requirements
# ------------

# Note: these should ultimately be better handled through a package structure
# etc.

# Required packages
library(distcrete)
library(epitrix)
library(here)
library(tidyverse) # total overkill, we merely use tibble and dplyr




# -------
# Helpers
# -------
# These are small functions to avoid repetitive code in the simulation

# Draw individual R based on group values and frequencies
# n: number of R values to draw
# draw_R <- function(n,
#                    group_R = c(R_detected, R_undetected),
#                    group_p = c(p_detected, 1 - p_detected)
# ) {
#   sample(
#     group_R,
#     size = n,
#     prob = group_p,
#     replace = TRUE)
# }


#------------------------------------------------
# Define the function. Its inputs are as follows:
#------------------------------------------------
#
# 1. R_undetected, the basic reproduction number when a case is undetected
#
# 2. intervention_efficacy, the prop reduction of R when cases detected 
#
# 3. serial_int, a distcrete object representing the serial interval
#    Use the 'make_disc_gamma' function to create this object
#
# 4. r_daily_intro, the daily rate of intro from reservoir
#
# 5. max_duration, the maximum duration of the outbreak
#
# 6. p_detected, proportion of cases detected


branching_process_model_mass <- function(  R_basic = 1.2, 
                                           intervention_efficacy = 0.9,
                                           serial_int_mean = 9, 
                                           serial_int_sd = 4,
                                           r_daily_intro = 50 / 784, 
                                           max_duration = 784,  
                                           n_sim = 1,
                                           vacc_coverage = 0.75,
                                           max_vac_eff = 0.9,
                                           vacc_time_to_max = 7,
                                           mean_delay_vacc2inf = 20,
                                           sd_delay_vacc2inf = 5,
                                           cluster_period = 30,
                                           reporting = 1.0,
                                           break_point = 365
                                         
) {
  
  source("make_disc_gamma.R")
  source("helper_functions_marburg_mass.R")
  source("create_logistic_function.R")
  
  serial_int <- make_disc_gamma(SI_mean, SI_sd)
  
  ## repro number when case is detected
  # R_detected <- R_undetected * (1 - intervention_efficacy) #<- without ring vaccination
  
  #R_detected 
  R_followed <- R_basic * (1 - intervention_efficacy)
  
  
  ## Create a vector of times for the loop
  vec_time <- seq_len(max_duration)
  
  # df_out <- data.frame(matrix(ncol=4,nrow=0, dimnames=list(NULL, c("n_sim", "cases", "introductions","max_ons_date"))))
  
  # --------------
  # Initialization
  # --------------
  
  # This is simulation stuff, but only done once at the beginning
  
  ## Introductions: we impose the first one on day 1, and others are drawn
  ## randomly
  
  n_daily_intro <- rpois(max_duration, r_daily_intro)
  #Ensure there is one introduced case on the first day of the outbreak
  n_daily_intro[1] <- 1L
  intro_onset <- rep(vec_time, n_daily_intro)
  
  library(tidyverse)
  ## Build a tibble - not absolutely needed, but nicer to inspect results
  # All cases will be stored in 'out'
  out <- tibble::tibble(case_id = seq_along(intro_onset),
                        date_onset = intro_onset)
  
  ## Determine R for each case
  
  # #Four types of cases
  # 
  # # 1: Naive case, not detected
  # R_undetected_unvacc <-  R_basic
  # 
  # # 2: Naive case, but detected
  # R_detected_unvacc <- R_basic * (1 - intervention_efficacy)
  # 
  # # 3: Vaccinated case, not detected
  # R_undetected_vacc <- R_basic * (1 - vacc_eff)
  # 
  # #4: Vaccinated, detected case
  # R_detected_vacc <- R_basic * (1 - intervention_efficacy) * (1 - vacc_eff)
  
  
  # Now add 2 more columns to 'out' - one specifying whether case happens 
  # before or after breakpoint (in words), the other translating this into a  
  # Boolean.
  
  out <- mutate(out,
                
                
                Timing = case_when(out$date_onset <= break_point ~ "Before breakpoint",
                                   out$date_onset > break_point ~ "After breakpoint"),
                
                p_detected = case_when(out$date_onset <= break_point ~ 0,
                                       out$date_onset > break_point ~ 1),
                
                
                
  ) 
  
  # Initialise 2 vectors - 1 that will store vaccine efficacies for each case
  # another to store the resultant reproduction number, R
  
  vac_eff <- vector()
  R_drawn <- vector()
  
  for (i in 1:length(out$case_id)) {
    
    vac_eff[i] =  make_logistic(max_vac_eff,
                                vacc_time_to_max,
                                mean_delay_vacc2inf,
                                sd_delay_vacc2inf,
                                cluster_period,
                                n_draws = 1)
    
    R_drawn[i] <- draw_R_mass(nrow(out),R_basic, intervention_efficacy,
                              vac_eff[i], out$p_detected[i],vacc_coverage)
  }
  
  # Add these to the new_cases dataframe
  out<- mutate(out, vac_eff=vac_eff, R = R_drawn)
  
  # Record the number of introductions: this is the number of rows in 'out'
  n_intros <- nrow(out)
  
  
  # -------------------------------------
  # Time iteration: what happens each day
  # -------------------------------------
  
  # Here we generate new cases from a branching process using the state of 'out'
  # above. The algorithm is:
  
  # 1. Determine the force of infection, i.e. the average number of new cases
  # generated.
  #
  # 2. Draw the number of new cases.
  #
  # 3. Draw features of the new cases.
  #
  # 4. Append new cases to the linelist of cases.
  #
  #
  # The main issue lies in determining the force of infection. For a single
  # individual 'i', it is determined as:
  #
  # lambda_i = R_i * w(t - onset_i)
  #
  # where 'w' is the PMF of the serial interval, 't' the current time, and
  # 'onset_i' the date of onset of case 'i'.
  #
  # The global force of infection at a given point in time is obtained by summing
  # all individual forces of infection.
  #
  
  
  # Time iteration
  for (t in 2:max_duration) {
    
    
    
    # Step 1
    lambda_i <- out$R * serial_int$d(t - out$date_onset)
    force_infection <- sum(lambda_i)
    
    # Step 2
    n_new_cases <- rpois(1, lambda = force_infection)
    
    # Step 3
    
    if(n_new_cases>0){
    last_id <- max(out$case_id)
    new_cases <- tibble(case_id = seq(from = last_id + 1,
                                      length.out = n_new_cases,
                                      by = 1L),
                        date_onset = rep(t, n_new_cases))

    
 
    # Now add 2 more columns to new_cases - one specifying whether case happens 
    # before or after breakpoint (in words), the other translating this into a  
    # Boolean.
    
    new_cases  <- mutate(new_cases ,
                  
                  
                  Timing = case_when(  new_cases$date_onset <= break_point ~ "Before breakpoint",
                                       new_cases$date_onset > break_point ~ "After breakpoint"),
                  
                  p_detected = case_when(  new_cases$date_onset <= break_point ~ 0,
                                           new_cases$date_onset > break_point ~ 1),
                  
    ) 
    
    # Initialise 2 vectors - 1 that will store vaccine efficacies for each case
    # another to store the resultant reproduction number, R
    
    vac_eff <- vector()
    R_drawn <- vector()
    
    for (i in 1:length(  new_cases$case_id)) {
      
      #Draw from the logistic curve
      vac_eff[i] =  make_logistic(max_vac_eff ,
                                  vacc_time_to_max,
                                  mean_delay_vacc2inf ,
                                  sd_delay_vacc2inf ,
                                  cluster_period ,
                                  n_draws = 1)
      
      # Use the helper function to draw a value of R, depending on vaccine params
      
      R_drawn[i] <- draw_R_mass(nrow(  new_cases ),R_basic, intervention_efficacy,
                                vac_eff[i],   new_cases$p_detected[i],vacc_coverage)
    }
    
    
    # Add these to the new_cases dataframe
    new_cases <- mutate(new_cases, vac_eff=vac_eff, R = R_drawn)
    
    
    # Step 4
    
    #Now add each new case to the 'out' dataframe (which stores all cases)
    out <- bind_rows(out, new_cases)
    
    }
    
    # Insert a break if we already have 5k cases
    if (nrow(out)>5000) {
      break
    }
  }
  
  # ------------
  # Check output
  # ------------
  
  # arrange cases by chronological order
  out <- arrange(out, date_onset)
  
  # redefine IDs to match chrono order
  out <- mutate(out, case_id = seq_len(n()))
  out

  
} #end function here


