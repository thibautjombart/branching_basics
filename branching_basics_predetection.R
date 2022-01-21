
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
source("create_logistic_function2.R")
source("helper_functions_marburg.R")



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


branching_process_with_detection <- function(
  R_basic = 1.2, 
  intervention_efficacy = 0.5,
  serial_int_mean = 9, 
  serial_int_sd = 4,
  r_daily_intro = 50 / 784, 
  max_duration = 784,  
  n_sim = 1,
  reporting = 0.8,
  break_point = 365
)
{
  
  # Build serial interval distribution
  serial_int <- make_disc_gamma(serial_int_mean, serial_int_sd)
  
  # R for detected cases
  
  
  # Define the different reproduction numbers
  
  # R_basic: R in absence of intervention
  # R_followed: R for cases who were followed (i.e. all cases but introductions)
  
  R_followed <- R_basic * (1 - intervention_efficacy)
  
  
  # Function to draw status of new cases
  
  # Note: this excludes introductions, which are determined at the beginning
  draw_status <- function(n = 1) {
    # Define probabilities of the different status
    
    # 'unreported': cases not reported
    # 'followed': cases reported and followed (i.e. non-introductions) 
    
    status_values <- c("unreported", "followed")
    status_prob <- c(1 - reporting, # unreported
                     reporting ) # followed
    
    sample(x = status_values,
           prob = status_prob,
           size = n,
           replace = TRUE)
  }
  
  
  # Create a vector of times for the loop
  vec_time <- seq_len(max_duration)
  
  
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
  
  ## Build a tibble - not absolutely needed, but nicer to inspect results
  out <- tibble::tibble(case_id = seq_along(intro_onset),
                        date_onset = intro_onset)
  
  ## Introductions are not impacted by follow-up so have max R if they occur pre-detection
  #Introduction after outbreak detection, though, have reduced R 
  out <- mutate(out,
                
                R = case_when( out$date_onset <= break_point ~ R_basic,
                               out$date_onset > break_point ~ R_followed),
                status = case_when( out$date_onset <= break_point ~ "intro, not followed",
                                    out$date_onset > break_point ~ "intro, followed")
                ) 
  
  
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
    
    
    # Technically only need to calculate R for cases on and before the timestep
    
    
    
    # Step 1
    lambda_i <- out$R * serial_int$d(t - out$date_onset)
    force_infection <- sum(lambda_i)
    
    # Step 2
    n_new_cases <- rpois(1, lambda = force_infection)
    
    # Step 3
    last_id <- max(out$case_id)
    
    new_cases <- tibble(case_id = seq(from = last_id + 1L,
                                      length.out = n_new_cases,
                                      by = 1L),
                        date_onset = rep(t, n_new_cases))
    
    # determine status: unreported or followed
    new_cases <- mutate(new_cases,
                        status = case_when(
                          new_cases$date_onset <= break_point ~ "not followed",
                          new_cases$date_onset > break_point ~ "followed"
                                           ),
                        R = case_when(
                            new_cases$date_onset <= break_point ~ R_basic,
                            new_cases$date_onset > break_point ~ R_followed
                                      )
                        )
    
    # Step 4
    out <- bind_rows(out, new_cases)
    
    
    # Insert a break if we already have >500 cases
    if (nrow(out) > 500) {
      break
    }
  }
  
  # arrange cases by chronological order
  out <- arrange(out, date_onset)
  
  # redefine IDs to match chrono order
  out <- mutate(out, case_id = seq_len(n()))
  out
} #end function here
