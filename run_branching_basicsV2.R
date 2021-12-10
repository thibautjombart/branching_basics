

# Required packages
library(distcrete)
library(epitrix)
library(here)
library(tidyverse) # total overkill, we merely use tibble and dplyr
library(foreach)
library(doParallel)
library(future.apply)
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
  wdir<-"C:/Users/grain/Dropbox/Marburg/Marburg/DRC_Ebola_HC_transmission_model-master/branching_basics"
  saved_files<-"C:/Users/grain/Dropbox/Marburg/Marburg/DRC_Ebola_HC_transmission_model-master"#"~/Documents/VEEPED/Ebola/Branching model with FOSA vac"
  setwd(wdir)
  outputs<-paste(saved_files,"/output_results_marburg/", sep="")
  # source(paste(wdir,"/ring_vaccine_multi_cluster_model_marburg.r", sep=""),local=TRUE) # Load model
  #  chains<-paste(saved_files,"/output_chain_marburg/",sep="")
}


SI_mean <- 9.2
SI_sd <- 4.4


p_detected <- 0.8
#Strategy: Try lower number of simulations in return for more parameter pairs
n_sim <- 50


# Outbreak duration and introduction rate are outbreak-specific parameters
# 
# #Frankfurt
# max_duration <- 100 # 784
# mean_intro_rate <- 20/max_duration
# max_cases <- 30
# min_cases<-10
# 
# 
# #Marburg
# max_duration <- 25
# mean_intro_rate <- 4/max_duration  
# max_cases <- 10
# min_cases<-0
# 
# 
# #Belgrade
# max_duration <- 17
# mean_intro_rate <- 1/max_duration 
# max_cases <- 10
# min_cases<-0
# 
# #SA
# max_duration <- 30
# mean_intro_rate <- 1/max_duration  
# max_cases <- 10
# min_cases<-0
# 
# 
# #Kenya
# max_duration <- 38
# mean_intro_rate <- 1/(max_duration+26) 
# max_cases <- 10
# min_cases<-0
# 
# 
# #Kenya2
# max_duration <- 10
# mean_intro_rate <- 1/(max_duration+26)  
# max_cases <- 10
# min_cases<-0
# 
# #DRC
max_duration <- 784
mean_intro_rate <- 50/max_duration
max_cases <- 200
min_cases<-100
# 
# #Angola
# max_duration <- 270
# mean_intro_rate <- 1/(max_duration+30) 
# max_cases <- 500
# min_cases<-200
# p_detected <- 0.5 # 1/3 of cases occurred before WHO intervention
# 
# #Uganda07
# max_duration <- 4
# mean_intro_rate <- 2/max_duration
# max_cases <- 10
# min_cases<-0
# p_detected <- 0.8
# 
# #Uganda08/09
# max_duration <- 180
# mean_intro_rate <- 2/max_duration
# max_cases <- 10
# min_cases<-0
# 
# #Uganda12
# max_duration <- 115
# mean_intro_rate <- 1/max_duration
# max_cases <- 50
# min_cases <- 10
# p_detected <- 0.3
# 
# #Uganda14
# max_duration <- 7
# mean_intro_rate <- 1/max_duration
# max_cases <- 10
# min_cases<-0
# p_detected <- 0.8
# 
# #Uganda17
# max_duration <- 44
# mean_intro_rate <- 1/max_duration
# max_cases <- 10
# min_cases<-0
# 


# Sample randomly from a grid

library(randtoolbox)

#Left column will be sample points for R, right column the sample points for int_efficacy
#UNCOMMENT THE FOLLWING 3 LINES IN FINAL VERSION
#sobol_draws <- sobol(n=3000, dim = 2, scrambling = 1 ,seed = 1880)

#write.csv(sobol_draws, "sobol_draws.csv", row.names=FALSE, quote=FALSE) 
#sobol_draws_save<- read.csv("sobol_draws.csv")

min_R_samp <- 0.50
max_R_samp <- 3.00
npoints<-1800

R_samp_location<-runif(n=npoints,min=min_R_samp,max=max_R_samp)
#R_samp_location <- sobol_draws_save[,1]*(max_R_samp-min_R_samp)+min_R_samp

points_sampled <- data.frame(R_und=R_samp_location,R_eff=runif(n=npoints))



#Chop off first part of points_sampled if necessary
# 
# index_counter <- 1800 #FILL THIS IN
# points_sampled<-points_sampled[-c(1:index_counter),]

#R_undetected_all <- seq(0.8,1.6, by=0.1)
#intervention_efficacy_all<- seq(0.6,1,by=0.1)

#draw from gamma distribution
r_daily_intro_all <- mean_intro_rate



#R_Int_Intro_grid<- expand.grid(R_undetected = R_undetected_all, intervention_efficacy = intervention_efficacy_all, r_daily_intro=r_daily_intro_all)

source("branching_basics2.R")


n_sim=1


#R_Int_Intro_df<-as.data.frame(R_Int_Intro_grid)

output <- future_mapply(branching_process_model,
              R_undetected=points_sampled[,1], # R_Int_Intro_df$R_undetected,
              intervention_efficacy=points_sampled[,2], #R_Int_Intro_df$intervention_efficacy, 
              serial_int_mean=SI_mean, 
              serial_int_sd=SI_sd,
              r_daily_intro=r_daily_intro_all,
              max_duration=max_duration,
              p_detected=p_detected,
              n_sim=n_sim)



df_output <- as.data.frame(output) 
n_cases_sim <- vector()

for (i in 1:length(df_output)) {
  
  n_cases_sim[i] <- df_output[i][3,]
  
  
}


df_results <- data.frame(points_sampled$R_und,points_sampled$R_eff ,n_cases_sim)

actual_cases <- 152

for (i in 1:length(df_output)){

    
    if( abs(df_results[i,3] - actual_cases) > 1/5*actual_cases )
      df_results[i,3] <- NA
  
}

Na_Status <- vector()

for (i in 1:length(df_output)){
  
  Na_Status[i] <- is.na(df_results[i,3]) 
  
}


Non_Na_Elems<- 1-Na_Status



R_Int_NonNAElems <- data.frame(df_results$points_sampled.R_und, df_results$points_sampled.R_eff,Non_Na_Elems) 

library(tidyverse)
library(MASS)
ggplot(R_Int_NonNAElems, aes(x=df_results$points_sampled.R_und, y=df_results$points_sampled.R_eff)) + 
  geom_point(aes(colour=factor(Non_Na_Elems)))+xlab("R for Undetected Cases") +ylab("Intervention Efficacy")


R_only_ones <- R_Int_NonNAElems %>% filter(Non_Na_Elems==1)

ggplot(R_only_ones, aes(x=df_results.points_sampled.R_und, y=df_results.points_sampled.R_eff)) + 
  geom_density_2d_filled(alpha=0.5)+xlab("R for Undetected Cases") +ylab("Intervention Efficacy")




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



