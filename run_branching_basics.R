

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
 
 
 p_detected <- 0.8
 n_sim <- 500
 
 
 # Outbreak duration and introduction rate are outbreak-specific parameters
 
 #Frankfurt
 max_duration <- 100 # 784
 mean_intro_rate <- 20/max_duration
 max_cases <- 30
 min_cases<-10
 
 
 #Marburg
 max_duration <- 25
 mean_intro_rate <- 4/max_duration  
 max_cases <- 10
 min_cases<-0
 
 
 #Belgrade
 max_duration <- 17
 mean_intro_rate <- 1/max_duration 
 max_cases <- 10
 min_cases<-0
 
 #SA
 max_duration <- 30
 mean_intro_rate <- 1/max_duration  
 max_cases <- 10
 min_cases<-0
 
 
 #Kenya
 max_duration <- 38
 mean_intro_rate <- 1/(max_duration+26) 
 max_cases <- 10
 min_cases<-0
 
 
 #Kenya2
 max_duration <- 10
 mean_intro_rate <- 1/(max_duration+26)  
 max_cases <- 10
 min_cases<-0
 
 #DRC
 max_duration <- 784
 mean_intro_rate <- 50/max_duration 
 max_cases <- 200
 min_cases<-100
 
 #Angola
 max_duration <- 270
 mean_intro_rate <- 1/(max_duration+30) 
 max_cases <- 500
 min_cases<-200
 
 #Uganda07
 max_duration <- 4
 mean_intro_rate <- 2/max_duration
 max_cases <- 10
 min_cases<-0
 
 #Uganda08/09
 max_duration <- 180
 mean_intro_rate <- 2/max_duration
 max_cases <- 10
 min_cases<-0
 
 #Uganda12
 max_duration <- 115
 mean_intro_rate <- 1/max_duration
 max_cases <- 50
 min_cases <- 10
 
 #Uganda14
 max_duration <- 7
 mean_intro_rate <- 1/max_duration
 max_cases <- 10
 min_cases<-0
 
 #Uganda17
 max_duration <- 44
 mean_intro_rate <- 1/max_duration
 max_cases <- 10
 min_cases<-0
 
 R_undetected_all <- seq(1.0,1.6, by=0.1)
 intervention_efficacy_all<- seq(0,1,by=0.2)
 
 #draw from gamma distribution
 r_daily_intro_all <- mean_intro_rate
 

 
 R_Int_Intro_grid<- expand.grid(R_undetected = R_undetected_all, intervention_efficacy = intervention_efficacy_all, r_daily_intro=r_daily_intro_all)
 
 source("branching_basics2.R")
 
 
 
 R_Int_Intro_df<-as.data.frame(R_Int_Intro_grid)
 
 future_mapply(branching_process_model,
               R_undetected=R_Int_Intro_df$R_undetected,
               intervention_efficacy=R_Int_Intro_df$intervention_efficacy, 
               r_daily_intro=R_Int_Intro_df$r_daily_intro,
               serial_int_mean=SI_mean, 
               serial_int_sd=SI_sd,
               max_duration)
 
 
 
#####################################################
 #Historical stuff
       
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



 