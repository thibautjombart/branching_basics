

# Required packages
library(distcrete)
library(epitrix)
library(here)
library(tidyverse) # total overkill, we merely use tibble and dplyr
library(foreach)
library(doParallel)
registerDoParallel(cores=3)


## WORKAROUND: https://github.com/rstudio/rstudio/issues/6692
## Revert to 'sequential' setup of PSOCK cluster in RStudio Console on macOS and R 4.0.0
if (Sys.getenv("RSTUDIO") == "1" && !nzchar(Sys.getenv("RSTUDIO_TERM")) && 
    Sys.info()["sysname"] == "Darwin" && getRversion() >= "4.0.0") {
  parallel:::setDefaultClusterOptions(setup_strategy = "sequential")
}

# Source external script to make discretized Gammas which we can use for PMF and
# for RNG
# path_to_file <- here("make_disc_gamma.R")
# source(path_to_file)
source("make_disc_gamma.R")
source("branching_basics2.R")

 
if(Sys.info()["user"]=="grain") {
  wdir<-"C:/Users/grain/Dropbox/Marburg/Marburg/DRC_Ebola_HC_transmission_model-master"
  saved_files<-"C:/Users/grain/Dropbox/Marburg/Marburg/DRC_Ebola_HC_transmission_model-master"#"~/Documents/VEEPED/Ebola/Branching model with FOSA vac"
  setwd(wdir)
  outputs<-paste(saved_files,"/output_results_marburg/", sep="")
  source(paste(wdir,"/ring_vaccine_multi_cluster_model_marburg.r", sep=""),local=TRUE) # Load model
  chains<-paste(saved_files,"/output_chain_marburg/",sep="")
}


 SI_mean <- 9.2
 SI_sd <- 4.4
 
 serial_int <- make_disc_gamma(SI_mean, SI_sd)


 max_duration <- 500
 
 p_detected <- 0.8
 n_sim <- 500
 
 R_undetected_all <- seq(1.2,2.0, by=0.1)
 R_detected_all <- seq(0.5,1.0, by=0.05)
 r_daily_intro_all <- seq(0,0.05, by=0.01)
 
 #foreach(i = c(1:length(R_undetected_all)),j=c(1:length(R_detected_all)),k=c(1:length(r_daily_intro_all)),.packages=c("tidyverse","doParallel","foreach","distcrete","epitrix") ) %dopar% { # Parallel code
  for (i in 1:length(R_undetected_all)) {
    for (j in 1:length(R_detected_all )) {
      for (k in 1:length(r_daily_intro_all)) {
       
       R_undetected <- R_undetected_all[i]
       R_detected <- R_detected_all[j]
       r_daily_intro <- r_daily_intro_all[k]
       
       
 intervention_efficacy <- (R_undetected-R_detected)/R_undetected
 
 
 source("branching_basics2.R")
 R_estimates <- branching_process_model(R_undetected, 
                                        intervention_efficacy, 
                                        R_detected,
                                        serial_int, 
                                        r_daily_intro, 
                                        max_duration, 
                                        p_detected,
                                        n_sim)

   }
  }
 }
 
# }

 