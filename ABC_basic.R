#ABC  

#Load relevant parameters for each outbreak

n_sim <- 500

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


#Load relevant results
setwd("C:/Users/grain/Dropbox/Marburg/Marburg/DRC_Ebola_HC_transmission_model-master/output_results_marburg/Angola")
temp_Angola = list.files(pattern="*.csv")
for (i in 1:length(temp_Angola)) assign(temp_Angola[i], read.csv(temp_Angola[i]))

n_cases_Angola <- data.frame(row.names=1:500)
n_intro_Angola <- data.frame(row.names=1:500)
for (i in 1:length(temp_Angola)){
  temp_df <- read.csv(temp_Angola[i])
  n_cases_Angola <- cbind(temp_df$V2, n_cases_Angola)
  colnames(n_cases_Angola)[1]<-temp_Angola[i]
  n_intro_Angola <- cbind(temp_df$V3, n_intro_Angola)
  colnames(n_intro_Angola)[1]<-temp_Angola[i]
}


# Put a 'NA' in place of predicted cases outside our tolerance region
#Angola: min: 200, max:500

for (i in 1:length(temp_Angola)){
  
  for (j in 1:n_sim){
    
    if( n_cases_Angola[i][j,] < min_cases | n_cases_Angola[i][j,] > max_cases)
    n_cases_Angola[i][j,] <- NA
    
    
    
  }
  
}


#Angola_Summary<- do.call(cbind, lapply(n_cases_Angola, summary))

#Now that we have NAs in the dataframe, we will need to use summarise_each
#and turn on the option: na.rm

Angola_Summary<- do.call(cbind, lapply(n_cases_Angola, summary,na.rm=TRUE))


#Also calculate how many predictions remain after we've filtered out
#unlikely estimates

Non_Na_Elems <-vector()

for (k in 1:length(temp_Angola)){

  #The 7th element of the summary stats shows the number of NAs
  Non_Na_Elems[k] <- n_sim - Angola_Summary[,k][7]  
  
}


# Plot a histogram of the posterior distribution

# Ignore NAs through !is.na

#n_cases_Angola_No_Na <- list()

for (k1 in 1:length(temp_Angola)){
  
  if(Non_Na_Elems[k1]==0){n_cases_Angola_No_Na<- -1}
  else{
n_cases_Angola_No_Na <- n_cases_Angola[k1][!is.na( n_cases_Angola[k1])]

hist(n_cases_Angola_No_Na,xlab="Number of Cases",ylab="Frequency",main= "Histogram of Cases ")

       }

}



#Do the above for other outbreaks 
