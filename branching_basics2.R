
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


branching_process_model <- function(R_undetected, 
                                    intervention_efficacy,
                                    R_detected,
                                    serial_int, 
                                    r_daily_intro, 
                                    max_duration, 
                                    p_detected,
                                    n_sim
) {
  
  ## repro number when case is detected
  # R_detected <- R_undetected * (1 - intervention_efficacy)
  
  ## Creat a vector of times for the loop
  vec_time <- seq_len(max_duration)
  
  df_out <- data.frame(matrix(ncol=4,nrow=0, dimnames=list(NULL, c("n_sim", "cases", "introductions","max_ons_date"))))
  
  # --------------
  # Initialization
  # --------------
  
  # This is simulation stuff, but only done once at the beginning
  
  ## Introductions: we impose the first one on day 1, and others are drawn
  ## randomly
  
  #for (count in 1:n_sim){
  
  library(foreach)
  library(doParallel)
  registerDoParallel(cores=2)
  
  ## WORKAROUND: https://github.com/rstudio/rstudio/issues/6692
  ## Revert to 'sequential' setup of PSOCK cluster in RStudio Console on macOS and R 4.0.0
  if (Sys.getenv("RSTUDIO") == "1" && !nzchar(Sys.getenv("RSTUDIO_TERM")) && 
      Sys.info()["sysname"] == "Darwin" && getRversion() >= "4.0.0") {
    parallel:::setDefaultClusterOptions(setup_strategy = "sequential")
  }
  
  df_result<-foreach(count = c(1:n_sim),.combine='cbind',.packages=c("tidyverse","doParallel")) %dopar% {
  
    
  n_daily_intro <- rpois(max_duration, r_daily_intro)
  n_daily_intro[1] <- 1L
  intro_onset <- rep(vec_time, n_daily_intro)
  
  library(tidyverse)
  ## Build a tibble - not absolutely needed, but nicer to inspect results
  out <- tibble::tibble(case_id = seq_along(intro_onset),
                date_onset = intro_onset)
  
  ## Determine R for each case
  source("helper_functions_marburg.R")
  out <- mutate(out,
                R = draw_R(nrow(out),R_detected,R_undetected,p_detected)
  )
  
  
  
  
  
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
    
    last_id <- max(out$case_id)
    new_cases <- tibble(case_id = seq(from = last_id + 1,
                                      length.out = n_new_cases,
                                      by = 1L),
                        date_onset = rep(t, n_new_cases))
    new_cases <- mutate(new_cases,
                        R = draw_R(n_new_cases,R_detected,R_undetected,p_detected)
    )
    
    # Step 4
    out <- bind_rows(out, new_cases)
  }
  
  # ------------
  # Check output
  # ------------
  library(incidence2)
  out %>%
    incidence(date_onset, interval = 7) %>%
    plot(title = "Simulated incidence", xlab = "Time")
  
  sim_number <- count
  n_cases <- nrow(out)
  n_intros <- as.data.frame(table(out$R))$Freq[1]
  max_onset_date <- max(out$date_onset)
  
  df_out[count, ] <- c(sim_number, n_cases ,n_intros,max_onset_date)
  
  }
  
  
  df_res_trans<-t(df_result)
  df_results<-as.data.frame(df_res_trans)
  
  df_results %>% 
    rename(
      sim_num = V1,
      num_cases = V2,
      num_intros = V3,
      maxim_ons_date = V4
    )
  
  #write.csv(out,"Results",R_undetected,round(R_detected,2),round(r_daily_intro,2),".csv",sep="")
  write.csv( df_results,file=paste0(outputs,"Results",
                                    "R_u=",R_undetected,"R_d=",round(R_detected,3),"Intro=",r_daily_intro*100,".csv"))
  
  # write.csv(output_df, file=paste0(outputs,"sims_ring_vac_only_P", 
  #                                  prop.ascertain*100, "_", paste("early_Guinea_Rmissed",format(round(r0_missed, 2), nsmall = 1),sep=""),
  #                                  paste("_Rwithin",format(round(r0_within_ring, 2), nsmall = 1),"_", sep=""),
  #                                  cap_max_days,"d_",n.sim,".csv"), row.names=F)
  # 
  #file=paste(chains,"cluster_chainsV",wvacc,"_P",100*pr.id,"_nrun",10000+nrun,".csv",sep="")
} #end function here


