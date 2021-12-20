
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
library(minpack.lm)

# Source external script to make discretized Gammas which we can use for PMF and
# for RNG
#path_to_file <- here("make_disc_gamma.R")
#source(path_to_file)
source("make_disc_gamma.R")



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


branching_process_model_ring <- function(
                                         R_undetected = 1.2, 
                                         intervention_efficacy = 0.75,
                                         serial_int_mean = 9, 
                                         serial_int_sd = 4,
                                         r_daily_intro = 50/784, 
                                         max_duration = 784,  
                                         p_detected = 0.8,
                                         n_sim = 1,
                                         cov_contacts = 0.5,
                                         max_vac_eff = 0.9 ,
                                         vacc_time_to_max = 12,
                                         vacc_delay_mean = 9.6,
                                         vacc_delay_sd = 4.2,
                                         serial_int_max = 21
                                
) {

  # Source external content
  source("make_disc_gamma.R")
  source("create_logistic_function2.R")
  source("helper_functions_marburg.R")
  
  serial_int <- make_disc_gamma(serial_int_mean, serial_int_sd)
  
  


  # R for detected cases
    
  # This version assumes a much longer period at max vaccine efficacy
  # Needs one more input - the max serial interval
  vac_eff <- make_logistic2(max_vac_eff,
                            vacc_time_to_max,
                            vacc_delay_mean,
                            vacc_delay_sd,
                            serial_int_max)
  
  
    ## repro number when case is detected
 
  R_detected                 <- ((1 - intervention_efficacy) *
                                 (1 - vac_eff*cov_contacts)) *
                                    R_undetected
  
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
  out <- tibble::tibble(case_id = seq_along(intro_onset),
                        date_onset = intro_onset)
  
  ## Determine R for each case

  out <- mutate(out,
                R = draw_R(nrow(out),R_detected,R_undetected,p_detected=0)
  )
  
  
  #Record the number of introductions: this is the number of rows in 'out'
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
    
    # Step 0: Modify R_detected
    vac_eff <- make_logistic2(max_vac_eff,
                              vacc_time_to_max,
                              vacc_delay_mean,
                              vacc_delay_sd,
                              serial_int_max)
    
    
    R_detected <- ((1 - intervention_efficacy) *
                     (1 - vac_eff*cov_contacts)) *
                      R_undetected
    
    #Technically only need to calculate R for cases on and before the timestep
    
    # out <- mutate(out,
    #               R = draw_R(nrow(out),R_detected,R_undetected,p_detected)
    # )
    # 
    
    # Step 1
    lambda_i <- out$R * serial_int$d(t - out$date_onset)
    force_infection <- sum(lambda_i)
    
    # Step 2
    n_new_cases <- rpois(1, lambda = force_infection)
    
    # Step 3
    
    last_id <- max(out$case_id)
    
    new_cases <- tibble(case_id = seq(from = last_id + 1,
                                      length.out = n_new_cases,
                                      by = 1L),
                        date_onset = rep(t, n_new_cases))
    
    new_cases <- mutate(new_cases,
                        R = draw_R(n_new_cases,R_detected,
                                   R_undetected,p_detected)
    )
    
    # Step 4
    out <- bind_rows(out, new_cases)
    
    
    # Insert a break if we already have >500 cases
    if (nrow(out)>500) {
      break
    }
  }
  
  # ------------
  # Check output
  # ------------
  
  
  # sim_number <- count
  
  #Total number of cases
  n_cases <- nrow(out)
  
  #Last day of outbreak
  max_onset_date <- max(out$date_onset)
  
  #df_out[count, ] <- c(sim_number, n_cases ,n_intros,max_onset_date)
  
  #library(incidence2)
  out
 } #end function here


