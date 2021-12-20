
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
path_to_file <- here("R", "make_disc_gamma.R")
source(path_to_file)



# ----------------------------------------------------------
# Parameters (hard coded, should be arguments of a function)
# ----------------------------------------------------------

## basic repro number
R_undetected <- 1.2

## prop reduction of R when cases detected 
intervention_efficacy <- 0.6

## repro number when case is detected
R_detected <- R_undetected * (1 - intervention_efficacy)

## reporting probability
reporting <- 0.8

## serial interval object
## this is a distcrete object
## $d is the pmf
## $r is the rng
serial_int <- make_disc_gamma(7, 3)

## Daily rate of intro from reservoir
r_daily_intro <- 0.05


## Maximum duration of the simulation
max_duration <- 365
vec_time <- seq_len(max_duration)

## Prop of cases detected
p_detected <- 0.7




# -------
# Helpers
# -------
# These are small functions to avoid repetitive code in the simulation

# Draw individual R based on group values and frequencies
# n: number of R values to draw

# note: here p_detected means 'of reported cases, which ones are attended in
# time for intervention to have efficacy'
draw_R <- function(n,
                   group_R = c(R_detected, R_undetected),
                   group_p = c(p_detected, 1 - p_detected)
                   ) {
  sample(
    group_R,
    size = n,
    prob = group_p,
    replace = TRUE)
}


# Draw status of cases: reported / unreported
draw_reported <- function(n = 1, p_reporting = 1) {
  sample(x = c("reported", "unreported"),
         prob = c(p_reporting, 1 - p_reporting),
         size = n,
         replace = TRUE)
}




# --------------
# Initialization
# --------------

# This is simulation stuff, but only done once at the beginning

## Introductions: we impose the first one on day 1, and others are drawn
## randomly

n_daily_intro <- rpois(max_duration, r_daily_intro)
n_daily_intro[1] <- 1L
intro_onset <- rep(vec_time, n_daily_intro)

## Build a tibble - not absolutely needed, but nicer to inspect results
out <- tibble(case_id = seq_along(intro_onset),
              date_onset = intro_onset)

## Determine R for each case
out <- mutate(out,
              R = draw_R(nrow(out)),
              status = draw_reported(nrow(out),
                                     p_reporting = reporting)
              )

## Make sure there is no intervention impact for unreported cases
out <- mutate(out,
              R = if_else(status == "unreported", R_undetected, R))



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
                      R = draw_R(n_new_cases),
                      status = draw_reported(nrow(new_cases),
                                             p_reporting = reporting)
                      )

  # Step 4
  out <- bind_rows(out, new_cases)

  ## Make sure there is no intervention impact for unreported cases
  out <- mutate(out,
                R = if_else(status == "unreported", R_undetected, R))

}




# ------------
# Check output
# ------------
library(incidence2)
out %>%
  incidence(date_onset, interval = 7) %>%
  plot(title = "Simulated incidence", xlab = "Time")
