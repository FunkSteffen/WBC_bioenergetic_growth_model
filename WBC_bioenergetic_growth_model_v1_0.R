###### WBC bioenergetic growth model
## Version: 1_0
## Author: Steffen Funk (2022)
## built in R version 4.1.2 (2021-11-01)
###

##set_working_directory ----
setwd("C:/Users/funk_/Desktop/WBC_Bioenergetic__model_v1_0/")

##read in required packages ----
library(ggplot2)
library(mgcv)
library(nnet)
library(MALDIquant)
library(plyr)
library(cowplot)


######LOAD ALL FUNCTIONS AND DATA NEEDED ----
#f1) temperature at residence depth ----
temperature_at_depth <- function(depth, j_day, y, temp_dat){
  allocated_temperature <- vector(length = length(depth), mode = "numeric")
  depth_strat_vector <- c(1.5,  4.5,  7.5, 10.5, 13.5, 16.5, 19.5, 22.5,
                          25.5, 28.5, 31.5, 34.5, 37.5, 40.5, 43.5, 46.5)
  for(i in 1:length(depth)){
    ind <- MALDIquant::match.closest(depth[i] , depth_strat_vector)
    chosen_strat <- depth_strat_vector[ind]
    sub_day <- subset(temp_dat, julian_day == j_day[i])
    sub_year <- subset(sub_day, year == y[i])
    sub_strat <- subset(sub_year, depth_strat == chosen_strat)
    allocated_temperature[i] <- sub_strat$mean_temp
  }
  return(allocated_temperature)
}
# with: depth - calculated residence depth, j_day -julian_day,
#       y - year, temp_dat - temperatures at 3m depth stratas taken from BSIOM(1979-2018)
#read in depending data set
temp_dat <- read.table("temp_at_3depth_strata_1979_2018.csv",header =TRUE, sep = ";")

#f2) SST and Tdiff for chosen model time series ----
temp_raw_dat <- function(y, raw_temp_dat){
  part1 <- subset(raw_temp_dat, year == y & julian_day  >= 46)
  part2 <- subset(raw_temp_dat, year == y+1 & julian_day  <= 45)
  full <- rbind(part1, part2)
  return(full)
}
# with: y - year, raw_temp_dat - SST and SBT data set taken from HH-Model (1979-2018)
# starts with 15.2. of chosen year and ends with 14.2. of year+1 (=> start with spent individuals and end up during spawning)
#read in depending data set
raw_temp_dat <- read.table("SST_SBT_Tdiff_1979_2018.csv",header =TRUE, sep = ";")

#f3) start length for modelling for a given age in a given year ----
length_start <- function(age, y, length_at_age_data){
  sub_dat <- subset(length_at_age_data, Year == y)
  sub_dat2 <- subset(sub_dat, Age == age)
  return(sub_dat2$mean_length)
}
# with age - chosen age, y - chosen age, y - chosen year, length_at_age_data -data set taken from bits smalk data

#f4) residence depth function ----
#load in depending rds data file including observer model
obs_model <- readRDS(file = "observer_model.rds")
residence_depth <- function(SST, Tdiff){
  new_dat <- data.frame(SST = SST, diff_sstsbt = Tdiff)
  pred <- predict(obs_model, newdata = new_dat, type = "response", se.fit = TRUE)
  pred_data <- data.frame(fit =pred$fit, random_resid = NA)
  for(i in 1:length(pred_data$fit)){
    pred_data$random_resid[i] <- sample(residuals(obs_model), size = 1, replace = TRUE)
  }
  residence_d_pred <- data.frame(residence_depth = pred_data$fit+ pred_data$random_resid)
  residence_d_pred$residence_depth[residence_d_pred$residence_depth < 0] <- 0
  residence_d_pred$residence_depth[residence_d_pred$residence_depth > 25] <- 25
  return(residence_d_pred)
}
#with: SST and tdiff chosen year (needs results from temp_raw_dat function)
#uses saved observer model

#f5) function to predict stomach content [g] ----
# calculates stomach content weight in g for chosen day for given length and at given residence depth
#load needed stomach content gam
scw_model <- readRDS(file = "stomach_model_incl_small.rds")
predict_stomach_content <- function(mean_depth, length, mean_temp){
  p <- predict(scw_model, newdata = data.frame(mean_depth = mean_depth,
                                                          length = length ,
                                                          mean_temp =mean_temp),
               se.fit = TRUE, type = "response")
  upr <- p$fit + (2 * p$se.fit)
  lwr <- p$fit - (2 * p$se.fit)
  scw_plus_random_effect <- p$fit + sample(residuals(scw_model),1,replace = TRUE)

  pred_scw <- exp(scw_plus_random_effect)
  
  ###correct with Mucus weight
  mucus_weight <- (1.692e-05 * length ^ 3.542e+00) * 0.1620566
  corrected_scw <- pred_scw - mucus_weight
  if(corrected_scw < 0 ){
    scw_final = 0
  } else {
    scw_final = corrected_scw
  }
  return(scw_final)
}

#example: calculate scw for first day of time series of chosen year 2016
#test_scw <-predict_stomach_content(mean_depth = depth_cod1_2016$residence_depth[1],
#                               length = 45, mean_temp = allocated_temperature_2016_cod1[1])
#plot(1:365,test_scw, type = "l")

#f6) function to define diet composition ----
diet_comp_mod <- readRDS(file = "multinom_mod.rds")
predict_dietclust <-function(month,depth, length){
  if(month >9){
    quarter <- 4
  }else{
      if(month > 6){
        quarter <- 3
      }else{
        if(month > 3){
          quarter <- 2
        }else{
      quarter <- 1
        }
      }
  }
  diet_clust_probs <- predict(diet_comp_mod,
                        newdata = data.frame(length = length,
                                             quarter = quarter,
                                             mean_depth =depth), type = "probs")
  diet_clust_array <- c(rep(1, times = round(diet_clust_probs[1]*100)),
                        rep(2, times = round(diet_clust_probs[2]*100)),
                        rep(3, times = round(diet_clust_probs[3]*100)),
                        rep(4, times = round(diet_clust_probs[4]*100)),
                        rep(5, times = round(diet_clust_probs[5]*100)),
                        rep(6, times = round(diet_clust_probs[6]*100)),
                        rep(7, times = round(diet_clust_probs[7]*100)),
                        rep(8, times = round(diet_clust_probs[8]*100)))
  diet_clust <- sample(diet_clust_array,1)
  return(diet_clust)
}
##test it:
#predicted_diet_cluster <- predict_dietclust(month = 1, depth = 20, length = 55)

#f7) daily consumption function ----
#calculates consumption per day in g
#read in pi_table

pi_table <- read.table("energy_dens_and_pi_table.csv", header = T, sep = ";") 

calc_consumption <- function(length, temp, scw, diet_clust){
  WF<- ps_lwr_a * length ^ ps_lwr_b
  psi <- pi_table$pi[pi_table$cluster_number == diet_clust]
  C24 <- 24*psi*(WF^0.305)*exp(0.11*temp)*scw^0.5
  return(C24)
}
##test it:
#cons <- calc_consumption(length = 45, temp = 15, scw = 5, diet_clust = 1)



#f8) daily maintenance ration function ----
calc_maint <- function(temp, diet_clust, activity_multiplier = 1.25, WF){
  prey_e_dens <- pi_table$energy_dens[pi_table$cluster_number == diet_clust]
  #WF<- ps_lwr_a * length ^ ps_lwr_b
  Rmaint <- (0.01203242* exp(0.056*temp)*(WF^0.736))* activity_multiplier/ (prey_e_dens*0.000239006) # prey e dens has to be calculated in Kcal
    return(Rmaint)
  }

##test it:
#maint <- calc_maint(temp =15, diet_clust = 1, activity_multiplier = 1.25, WF = 800)

#f9) calculate food for growth as function ----
#calculates food for growth in g per day
calc_ffg <-function(C24, rmaint){
  ffg <- C24-rmaint
  return(ffg)
}
##test it:
#ffg <-calc_ffg(17,10)

#f10) weight gain function ----
#uses ffg and overall net conversion efficiency estimate of 0.35 (Temming and Herrmann, 2009)
weight_gain <- function(ffg){
  gain <- ffg*0.35
  return(gain)
}


#f11) numeric growth as function ----
grow <- function(start_length, chosen_year){
  
  #define length weight relationship coefficients
  ps_lwr_a = 0.01052065
  ps_lwr_b = 2.969391
  s_lwr_a = 0.007808755
  s_lwr_b = 3.049126
  
  start_weight <-  ps_lwr_a *start_length ^ ps_lwr_b
  chosen_temp <- temp_raw_dat(y = chosen_year, raw_temp_dat)
  new_matrix <- matrix(data = NA, nrow = nrow(chosen_temp), ncol = 19)
  new_matrix[,1] <-chosen_temp$year
  new_matrix[,2] <- chosen_temp$julian_day
  new_matrix[,3] <- chosen_temp$month
  new_matrix[,4] <- chosen_temp$SST
  new_matrix[,5] <- chosen_temp$tdiff
  new_matrix[,6] <- unlist(residence_depth(SST = new_matrix[,4], new_matrix[,5]))
  new_matrix[,7] <- unlist(temperature_at_depth(depth = new_matrix[,6],
                                                j_day = new_matrix[,2],
                                                y = new_matrix[,1],
                                                temp_dat = temp_dat))
  new_matrix[1,8] <- start_length # start length
  new_matrix[1,9] <- start_weight # start weight
  new_matrix[1,10] <- ps_lwr_a *new_matrix[1,8] ^ ps_lwr_b # weight for physiologic processes -cons and Rmaint
  new_matrix[1,11] <- predict_stomach_content(mean_depth = new_matrix[1,6],
                                              length = new_matrix[1,8], mean_temp = new_matrix[1,7])
  new_matrix[1,12] <- predict_dietclust(month = new_matrix[1,3], depth = new_matrix[1,6],
                                        length = new_matrix[1,8])
  new_matrix[1,13] <- calc_consumption(length = new_matrix[1,8], temp = new_matrix[1,7],
                                       scw = new_matrix[1,11], diet_clust = new_matrix[1,12])
  #activity multiplier 
  if(new_matrix[1,13] > 0){
    act <- 1.25
  }else{
    act <- 1
  }
  new_matrix[1,14] <- calc_maint(temp = new_matrix[1,7], diet_clust = new_matrix[1,12],
                                 activity_multiplier = act, WF = new_matrix[1,9])
  new_matrix[1,15] <- calc_ffg(C24 = new_matrix[1,13], rmaint = new_matrix[1,14])
  new_matrix[1,16] <- weight_gain(new_matrix[1,15])
  new_matrix[1,17] <- new_matrix[1,9] + new_matrix[1,16]
  new_matrix[1,18] <- (new_matrix[1,17]/s_lwr_a)^(1/s_lwr_b) # new pred length
  new_matrix[1,19] <- 1
  
  for(i in 2:nrow(new_matrix)){
    #if start length is longer than pred length, start length will be used as length
    if(new_matrix[i-1,8] > new_matrix[i-1,18]){
      new_matrix[i,8] <- new_matrix[i-1,8]
    }else{
      new_matrix[i,8] <- new_matrix[i-1,18]
    }
    #new start weight is end weight of previous step
    new_matrix[i,9] <- new_matrix[i-1,17]
    #physiological weight is fixed to actual length of the cod
    new_matrix[i,10] <- ps_lwr_a *new_matrix[i,8] ^ ps_lwr_b
    #scw
    new_matrix[i,11] <- predict_stomach_content(mean_depth = new_matrix[i,6],
                                                length = new_matrix[i,8], mean_temp = new_matrix[i,7])
    #diet clust
    new_matrix[i,12] <- predict_dietclust(month = new_matrix[i,3], depth = new_matrix[i,6],
                                          length = new_matrix[i,8])
    #C24
    new_matrix[i,13] <- calc_consumption(length = new_matrix[i,8], temp = new_matrix[i,7],
                                         scw = new_matrix[i,11], diet_clust = new_matrix[i,12])
    #activity multiplier 
    if(new_matrix[i,13] > 0){
      act <- 1.25
    }else{
      act <- 1
    }
    #rmaint
    new_matrix[i,14] <- calc_maint(temp = new_matrix[i,7],
                                   diet_clust = new_matrix[i,12], activity_multiplier = act,
                                   WF = new_matrix[i,9]) # related to actual weight of the cod
    #ffg
    new_matrix[i,15] <- calc_ffg(C24 = new_matrix[i,13], rmaint = new_matrix[i,14])
    #weight gain
    new_matrix[i,16] <- weight_gain(new_matrix[i,15])
    #new weight
    new_matrix[i,17] <- new_matrix[i,9] + new_matrix[i,16]
    #new predicted length
    new_matrix[i,18] <- (new_matrix[i,17]/s_lwr_a)^(1/s_lwr_b)
    new_matrix[i,19] <- i
  }
  return(new_matrix)
}

#f12) grow and plotting as function ----
#smalk_data <- read.table("smalk_raw_data_males.csv",header = TRUE, sep = ";")
#number_ind = 10
#chosen_year = 2016
#chosen_age = 3

grow_and_plot <- function(chosen_year, chosen_age, number_ind, smalk_data = smalk_data){
  chosen_year <- chosen_year
  chosen_age <- chosen_age
  number_ind <- number_ind
  end_length_vector <- sample(smalk_data$length[smalk_data$year == chosen_year+1 & 
                                                  smalk_data$age == chosen_age+1]/10, number_ind, 
                              replace = TRUE)
  n_end_indiv <- length(smalk_data$length[smalk_data$year == chosen_year+1 & 
                                            smalk_data$age == chosen_age+1])
  start_length_vector <- smalk_data$length[smalk_data$year == chosen_year & smalk_data$age == chosen_age]/10
  n_start_indiv <- length(smalk_data$length[smalk_data$year == chosen_year & smalk_data$age == chosen_age])
  #use grow function
  
  new_list <- vector(mode = "list", length = number_ind)
  for(i in 1:number_ind){
    start_length <- sample(start_length_vector,1,replace = TRUE) # random start length
    test1 <- grow(start_length = start_length, chosen_year = chosen_year)
    new_list[[i]] <-test1
  }
  #bind all matrices
  j <- do.call("rbind",new_list)
  j_data <- as.data.frame(j)
  names(j_data) <- c("year","julian_day","month","SST","tdiff","res_depth","temp_at_depth","start_length",
                     "start_weight","phys_weight","scw","diet_clust","consumption","rmaint","ffg",
                     "weight_gain","new_weight",
                     "pred_length","growth_day")
  
  check2 <- ddply(j_data, .(growth_day), summarise, med_resdepth= median(res_depth), perc90 = quantile(res_depth,.90),
                  perc10 = quantile(res_depth,.10))
  residence_depth_plot2 <-ggplot(data = check2, aes(x = growth_day, y = med_resdepth))+geom_line()+
    geom_ribbon(aes(ymin=perc10, 
                    ymax=perc90), alpha=.2)+theme_bw()+
    ylab("Residence depth [m]") + xlab("modelling time [days]") + ylim(35,0)
  
  ##plot daily consumption in g
  j_data$rel_cons <- j_data$consumption/j_data$start_weight
  j_data$rel_maint <- j_data$rmaint/j_data$start_weight
  consum <- ddply(j_data, .(growth_day), summarise, med_cons= median(rel_cons), cons10 = quantile(rel_cons,.10),
                  cons90 = quantile(rel_cons,.90),
                  med_maint = median(rel_maint), 
                  maint90 = quantile(rel_maint,.90), maint10 = quantile(rel_maint,.10))
  
  rel_cons_plot <- ggplot(data = consum, aes(x = growth_day, y = med_cons * 100))+geom_line(col = "blue")+
    geom_ribbon(aes(ymin= cons10*100,
                    ymax= cons90*100), alpha=.2, fill = "blue")+theme_bw()+
    ylab("consumption [% body weight]") + xlab("modelling time [days]") 
  
  rel_maint_plot <- ggplot(data = consum, aes(x = growth_day, y = med_maint * 100))+geom_line(col = "red")+
    geom_ribbon(aes(ymin=maint10*100,
                    ymax=maint90*100), alpha=.2, fill = "red")+theme_bw()+
    ylab("Maintenance ration [% body weight]") + xlab("modelling time [days]") 
  
  ##plot ffg in percent body weight over time
  j_data$ffg_rel <- j_data$ffg/j_data$start_weight
  check_ffg <- ddply(j_data, .(growth_day), summarise, med_ffg= median(ffg_rel), ffg90 = quantile(ffg_rel,.90),
                     ffg10 = quantile(ffg_rel,.10))
  
  rel_ffg_plot <- ggplot(data = check_ffg, aes(x = growth_day, y = med_ffg * 100))+geom_line()+
    geom_ribbon(aes(ymin=ffg10*100,
                    ymax=ffg90*100), alpha=.2)+theme_bw()+
    ylab("food for growth [% body weight]") + xlab("modelling time [days]") 
  
  bodyweight <- ddply(j_data, .(growth_day), summarise, med_bw = median(start_weight),
                      weight10 = quantile(start_weight, .10), weight90 = quantile(start_weight, .90))
  bodyweight_plot <- ggplot(data = bodyweight, aes(x = growth_day, y = med_bw))+geom_line()+
    geom_ribbon(aes(ymin=weight10,
                    ymax=weight90), alpha=.2)+theme_bw()+
    ylab("body weight [g]") + xlab("modelling time [days]") 
  
  boxplot_end_length <- ggplot(data =rbind(data.frame(length = floor(j_data$pred_length[j_data$julian_day == 45]), 
                                                      group = "predicted"),
                                           data.frame(length = end_length_vector,  group = "observed")), 
                               aes(x = group, y = length))+geom_boxplot()+
    xlab("")+ylab("length [cm]") + theme_bw()
  
  pred_plot <- cowplot::plot_grid(residence_depth_plot2, rel_cons_plot , rel_maint_plot, 
                                  rel_ffg_plot, bodyweight_plot, boxplot_end_length, labels = "auto")
  
  return(pred_plot)
  
}
#test <- grow_and_plot(chosen_year = 2016, chosen_age = 3, number_ind = 25,
#                      smalk_data = smalk_data)

####EXAMPLES ----

###### example 1: predict residence depth of a cod in 2016 and plot temperature at residence depth----
#calculate tdiff and SSt for selcted year 2016
temp2016 <- temp_raw_dat(y = 2016, raw_temp_dat = raw_temp_dat)

#predict residence depth for a cod in 2016
depth_cod1_2016 <- residence_depth(SST = temp2016$SST, Tdiff = temp2016$tdiff)
#depth_cod1_2016b <- residence_depth2(SST = temp2016$SST, Tdiff = temp2016$tdiff)
plot(1:365, depth_cod1_2016$residence_depth, type = "l")
#lines(1:365, depth_cod1_2016b$residence_depth, col = 2)

#allocate temperature to residence depth
allocated_temperature_2016_cod1 <- temperature_at_depth(depth = depth_cod1_2016$residence_depth,
                                                        j_day = temp2016$julian_day,
                                                        y = temp2016$year,
                                                        temp_dat = temp_dat)

plot(1:365, allocated_temperature_2016_cod1)



###### example 2: numeric growth of a single male cod of age 3 in 2016-2017 ----

#read in length at age data set
smalk_data <- read.table("smalk_raw_data_males.csv",header = TRUE, sep = ";")

chosen_year <- 2016

end_length_vector <- rnorm(n = 25, mean = mean(smalk_data$length[smalk_data$year == chosen_year+1 & smalk_data$age == 4]),
                      sd = sd(smalk_data$length[smalk_data$year == chosen_year+1 & smalk_data$age == 4]))
start_length <- smalk_data$length[smalk_data$year == chosen_year & smalk_data$age == 3]

#define length weight relationship coefficients
ps_lwr_a = 0.01052065
ps_lwr_b = 2.969391
s_lwr_a = 0.007808755
s_lwr_b = 3.049126

start_weight <- ps_lwr_a *(start_length/10) ^ ps_lwr_b
phys_weight <- ps_lwr_a *(start_length/10) ^ ps_lwr_b

chosen_temp <- temp_raw_dat(y = chosen_year, raw_temp_dat)

#growth matrix
new_matrix <- matrix(data = NA, nrow = nrow(chosen_temp), ncol = 18)
new_matrix[,1] <-chosen_temp$year
new_matrix[,2] <- chosen_temp$julian_day
new_matrix[,3] <- chosen_temp$month
new_matrix[,4] <- chosen_temp$SST
new_matrix[,5] <- chosen_temp$tdiff
new_matrix[,6] <- unlist(residence_depth(SST = new_matrix[,4], new_matrix[,5]))
new_matrix[,7] <- unlist(temperature_at_depth(depth = new_matrix[,6],
                                       j_day = new_matrix[,2],
                                       y = new_matrix[,1],
                                       temp_dat = temp_dat))
new_matrix[1,8] <- start_length[1]/10 # start length
new_matrix[1,9] <- start_weight[1] # start weight
new_matrix[1,10] <- ps_lwr_a *new_matrix[1,8] ^ ps_lwr_b # weight for physiologic processes -cons and Rmaint
new_matrix[1,11] <- predict_stomach_content(mean_depth = new_matrix[1,6],
                                            length = new_matrix[1,8], mean_temp = new_matrix[1,7])
new_matrix[1,12] <- predict_dietclust(month = new_matrix[1,3], depth = new_matrix[1,6],
                                      length = new_matrix[1,8])
new_matrix[1,13] <- calc_consumption(length = new_matrix[1,8], temp = new_matrix[1,7],
                                     scw = new_matrix[1,11], diet_clust = new_matrix[1,12])
new_matrix[1,14] <- calc_maint(temp = new_matrix[1,7], diet_clust = new_matrix[1,12],
                               activity_multiplier = 1.25, WF = new_matrix[1,9])
new_matrix[1,15] <- calc_ffg(C24 = new_matrix[1,13], rmaint = new_matrix[1,14])
new_matrix[1,16] <- weight_gain(new_matrix[1,15])
new_matrix[1,17] <- new_matrix[1,9] + new_matrix[1,16]
new_matrix[1,18] <- (new_matrix[1,17]/s_lwr_a)^(1/s_lwr_b) # new pred length

for(i in 2:nrow(new_matrix)){
  #if start length is longer than pred length, start length will be used as length
  if(new_matrix[i-1,8] > new_matrix[i-1,18]){
    new_matrix[i,8] <- new_matrix[i-1,8]
  }else{
    new_matrix[i,8] <- new_matrix[i-1,18]
  }
  #new start weight is end weight of previous step
  new_matrix[i,9] <- new_matrix[i-1,17]
  #physiological weight is fixed to actual length of the cod
  new_matrix[i,10] <- ps_lwr_a *new_matrix[i,8] ^ ps_lwr_b
  #scw
  new_matrix[i,11] <- predict_stomach_content(mean_depth = new_matrix[i,6],
                                              length = new_matrix[i,8], mean_temp = new_matrix[i,7])
  #diet clust
  new_matrix[i,12] <- predict_dietclust(month = new_matrix[i,3], depth = new_matrix[i,6],
                                        length = new_matrix[i,8])
  #C24
  new_matrix[i,13] <- calc_consumption(length = new_matrix[i,8], temp = new_matrix[i,7],
                                       scw = new_matrix[i,11], diet_clust = new_matrix[i,12])
  #rmaint
  new_matrix[i,14] <- calc_maint(temp = new_matrix[i,7],
                                 diet_clust = new_matrix[i,12], activity_multiplier = 1.25,
                                 WF = new_matrix[i,9]) # auf aktuelles Gewicht bezogen
  #ffg
  new_matrix[i,15] <- calc_ffg(C24 = new_matrix[i,13], rmaint = new_matrix[i,14])
  #weight gain
  new_matrix[i,16] <- weight_gain(new_matrix[i,15])
  #new weight
  new_matrix[i,17] <- new_matrix[i,9] + new_matrix[i,16]
  #new predicted length
  new_matrix[i,18] <- (new_matrix[i,17]/s_lwr_a)^(1/s_lwr_b)
}

plot(1:365,new_matrix[,13], type = "l")
lines(1:365,new_matrix[,14], col =2)

#plot ffg over time
plot(new_matrix[,2],new_matrix[,15]/new_matrix[,17], type = "l")

#weight over time
plot(1:365, new_matrix[,17])

#end_length
new_matrix[nrow(new_matrix), 18]

#stomach content, daily consumption and daily maintenance ration over the modelling period
plot(1:365,new_matrix[,11], type = "l") # stomach content [g]
lines(1:365,new_matrix[,13], col = 2) # daily consumption [g]
lines(1:365,new_matrix[,14], col =3) # daily maintenance ration [g]



######example 3: usage of "grow"-function and subsequently plotting (for 10 male individuals in 2016) ----
smalk_data <- read.table("smalk_raw_data_males.csv",header = TRUE, sep = ";")

#define length weight relationship coefficients
ps_lwr_a = 0.01052065
ps_lwr_b = 2.969391
s_lwr_a = 0.007808755
s_lwr_b = 3.049126

#choose modelling year
chosen_year <- 2017

#choose which year class should be modelled
chosen_age <- 3

#export start and end length values (observed data for later comparison) from corresponding BITS data
end_length_vector <- sample(smalk_data$length[smalk_data$year == chosen_year+1 & smalk_data$age == chosen_age+1]/10, 100, replace = TRUE)
hist(end_length_vector)

start_length_vector <- smalk_data$length[smalk_data$year == chosen_year & smalk_data$age == chosen_age]/10
hist(start_length_vector)

new_list <- vector(mode = "list", length = 10)
for(i in 1:10){
  start_length <- sample(start_length_vector,1,replace = TRUE) # random start length
  test1 <- grow(start_length = start_length, chosen_year = chosen_year)
  new_list[[i]] <-test1
}

#prepare data for subsequently plotting
#bind all matrices
j <- do.call("rbind",new_list)

j_data <- as.data.frame(j)
names(j_data) <- c("year","julian_day","month","SST","tdiff","res_depth","temp_at_depth","start_length",
                   "start_weight","phys_weight","scw","diet_clust","consumption","rmaint","ffg","weight_gain","new_weight",
                   "pred_length","growth_day")
summary(j_data)

##calculate it with 5 and 95% percentile instead of standard deviation
check2 <- ddply(j_data, .(growth_day), summarise, med_resdepth= median(res_depth), perc90 = quantile(res_depth,.90),
                perc10 = quantile(res_depth,.10))
residence_depth_plot2 <-ggplot(data = check2, aes(x = growth_day, y = med_resdepth))+geom_line()+
  geom_ribbon(aes(ymin=perc10, 
                  ymax=perc90), alpha=.2)+theme_bw()+
  ylab("Residence depth [m]") + xlab("modelling time [days]") + ylim(35,0)

##plot daily consumption in g
j_data$rel_cons <- j_data$consumption/j_data$start_weight
j_data$rel_maint <- j_data$rmaint/j_data$start_weight
consum <- ddply(j_data, .(growth_day), summarise, med_cons= median(rel_cons), cons10 = quantile(rel_cons,.10),
                cons90 = quantile(rel_cons,.90),
                med_maint = median(rel_maint), 
                maint90 = quantile(rel_maint,.90), maint10 = quantile(rel_maint,.10))


rel_cons_plot <- ggplot(data = consum, aes(x = growth_day, y = med_cons * 100))+geom_line(col = "blue")+
  geom_ribbon(aes(ymin= cons10*100,
                  ymax= cons90*100), alpha=.2, fill = "blue")+theme_bw()+
  ylab("consumption [% body weight]") + xlab("modelling time [days]") 

rel_maint_plot <- ggplot(data = consum, aes(x = growth_day, y = med_maint * 100))+geom_line(col = "red")+
  geom_ribbon(aes(ymin=maint10*100,
                  ymax=maint90*100), alpha=.2, fill = "red")+theme_bw()+
  ylab("Maintenance ration [% body weight]") + xlab("modelling time [days]") 

##plot ffg in percent body weight over time
j_data$ffg_rel <- j_data$ffg/j_data$start_weight
check_ffg <- ddply(j_data, .(growth_day), summarise, med_ffg= median(ffg_rel), ffg90 = quantile(ffg_rel,.90),
                   ffg10 = quantile(ffg_rel,.10))

rel_ffg_plot <- ggplot(data = check_ffg, aes(x = growth_day, y = med_ffg * 100))+geom_line()+
  geom_ribbon(aes(ymin=ffg10*100,
                  ymax=ffg90*100), alpha=.2)+theme_bw()+
  ylab("food for growth [% body weight]") + xlab("modelling time [days]") 

boxplot_end_length <- ggplot(data =rbind(data.frame(length = j_data$pred_length[j_data$julian_day == 45], 
                                                    group = "predicted"),
                                         data.frame(length = end_length_vector,  group = "observed")), 
                             aes(x = group, y = length))+geom_boxplot()+
  xlab("")+ylab("length [cm]") + theme_bw()
print(boxplot_end_length)

plot_grid(residence_depth_plot2, rel_cons_plot , rel_maint_plot, 
                   rel_ffg_plot, boxplot_end_length, labels = "auto")


######example 4: run and directly plot models for males ----
smalk_data <- read.table("smalk_raw_data_males.csv",header = TRUE, sep = ";")
test <- grow_and_plot(chosen_year = 2016, chosen_age = 3, number_ind = 10,
                      smalk_data = smalk_data)

test #print plot_grid
