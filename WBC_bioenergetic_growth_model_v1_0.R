######
## growth model v.1.0
## setup model setting
## author: Steffen Funk (2019)
######

##set_working_directory ----
#setwd("~/")

###load in allocation functions ----
#f1) temperature at residence depth
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
#       y - year, temp_dat - temperatures at 3m depth stratas taken from HH-Model (1979-2018)
#read in depending data set
temp_dat <- read.table("temp_at_3depth_strata_1979_2018.csv",header =TRUE, sep = ";")

#f2) SST and tdiff for modelled time series
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

#f3) start length for modelling for a given age in a given year
length_start <- function(age, y, length_at_age_data){
  sub_dat <- subset(length_at_age_data, Year == y)
  sub_dat2 <- subset(sub_dat, Age == age)
  return(sub_dat2$mean_length)
}
# with age - chosen age, y - chosen age, y - chosen year, length_at_age_data -data set taken from bits smalk data
# Areas SD 22, 23,24, all maturities and length and weight calculated as mean (only females)
#read in depending data set
length_at_age_data <- read.table("length_weight_at_age_females_13_19.csv",header =TRUE, sep = ";")


##define coeffients for length weight_relation_ships ----

## length weight relationship coefficients calcullated with smalk data (only females) for post spawners (ps)
## over Sd22-23, 1991-2019 (R^2 = 0.963)
## use these coefficients to caluclate start weight (After spawning!)
ps_lwr_a <- 0.01361447
ps_lwr_b <- 2.91643490

## length weight relationship coefficients calcullated with smalk data (only females) for spawners (s)
## over Sd22-23, 1991-2019 (R^2 = 0.950)
s_lwr_a <- 0.006041833
s_lwr_b <- 3.144840058


#### calculate residence depth
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

#alternative use of residence depth prediction (not adding random residual, draw a value between confidence intervals instead)
#residence_depth2 <- function(SST, Tdiff){
#  new_dat <- data.frame(SST = SST, diff_sstsbt = Tdiff)
#  pred <- predict(obs_model, newdata = new_dat, type = "response", se.fit = TRUE)
#  resd_uppr <- pred$fit + pred$se.fit *1.96
#  resd_lwr <- pred$fit - pred$se.fit *1.96
#  depth_vec <- vector(mode = "list", length =length(pred$se.fit))
#  for(i in 1: length(pred$se.fit)){
#    depth_vec[[i]] <- sample(seq(from = resd_lwr[i], to = resd_uppr[i], length.out = 15),1,replace = TRUE)
#  }
#  residence_d_pred <- data.frame(residence_depth = unlist(depth_vec))
#  residence_d_pred$residence_depth[residence_d_pred$residence_depth < 0] <- 0
#  residence_d_pred$residence_depth[residence_d_pred$residence_depth > 25] <- 25
#  return(residence_d_pred)
#}




#with: SST and tdiff chosen year (needs results from temp_raw_dat function)
#uses saved observer model

####calculate stomach content [g]
# calculates stomach content weight in g for chosen day for given length and at given residence depth
library(mgcv)
#load needed stomach content gam
#scw_model <- readRDS(file = "stomach_model.rds")
scw_model <- readRDS(file = "stomach_model_incl_small.rds")
predict_stomach_content <- function(mean_depth, length, mean_temp){
  p <- predict(scw_model, newdata = data.frame(mean_depth = mean_depth,
                                                          length = length ,
                                                          mean_temp =mean_temp),
               se.fit = TRUE, type = "response")
  upr <- p$fit + (2 * p$se.fit)
  lwr <- p$fit - (2 * p$se.fit)
  #scw_plus_random_effect <- p$fit + sample(residuals(scw_model),1,replace = TRUE)
  
  scw_plus_random_effect <- sample(seq(lwr, upr, length.out = 10),1, replace = TRUE)

  pred_scw <- exp(scw_plus_random_effect)
  
  ###correct with Mucus weight
  ####korrigiere stomach content weight mit Mucus
  # Mucus-Menge bleibt konstant über die Fischlänge und Mucus Masse [g] 
  # entspricht immer ungefähr  data_Mucus$stomach_empty[g]* 0.1620566
  # Magenwandgewicht [g]= 1.692e-05 * length [cm] ^ 3.542e+00 
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
test_scw <-predict_stomach_content(mean_depth = depth_cod1_2016$residence_depth[1],
                               length = 45, mean_temp = allocated_temperature_2016_cod1[1])
#plot(1:365,test_scw, type = "l")

####calcualte diet composition
library(nnet)
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
  diet_clust <- predict(diet_comp_mod,
                        newdata = data.frame(length = length,
                                             quarter = quarter,
                                             mean_depth =depth))
  return(diet_clust)
}
##test it:
#predicted_diet_cluster <- predict_dietclust(month = 1, depth = 20, length = 55)

####calculate consumption
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

####calculate Maintenance ration
#calculates maintenance ration per day in g (taking espect to the diet composition by using comp-specific energy dens)
#calc_maint <- function(length, temp, diet_clust, activity_multiplier = 1, WF){
#  #WF<- ps_lwr_a * length ^ ps_lwr_b
#  prey_e_dens <- pi_table$energy_dens[pi_table$cluster_number == diet_clust]
#  Rmaint <- (0.0097* exp(0.056*temp)*(WF^0.76))* activity_multiplier/ (prey_e_dens*0.000239006) # prey e dens has to be calculated in Kcal
#  return(Rmaint)
#}

###calc maintenance via Panten and Saunders
calc_maint <- function(temp, diet_clust, activity_multiplier = 1.25, WF){
  prey_e_dens <- pi_table$energy_dens[pi_table$cluster_number == diet_clust]
  #WF<- ps_lwr_a * length ^ ps_lwr_b
  Rmaint <- (0.01203242* exp(0.056*temp)*(WF^0.736))* activity_multiplier/ (prey_e_dens*0.000239006) # prey e dens has to be calculated in Kcal
    return(Rmaint)
  }

##test it:
#maint <- calc_maint(temp =15, diet_clust = 1, activity_multiplier = 1.25, WF = 800)

####calculate food for growth
#calculates food for growth in g per day
calc_ffg <-function(C24, rmaint){
  ffg <- C24-rmaint
  return(ffg)
}
##test it:
#ffg <-calc_ffg(17,10)

####calculate weight gain
#uses ffg and overall net conversion efficiency estimate of 0.55 (Temming and Herrmann, 2009)
weight_gain <- function(ffg){
  gain <- ffg*0.55
  return(gain)
}

#### numeric growth as function ----
grow <- function(start_length, chosen_year){
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
                                   WF = new_matrix[i,9]) # auf aktuelles Gewicht bezogen
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
# let one fish grow individually (uses all functions above)

#grow2 <- function(start_length, chosen_year){
  start_weight <-  ps_lwr_a *start_length ^ ps_lwr_b
  chosen_temp <- temp_raw_dat(y = chosen_year, raw_temp_dat)
  new_matrix <- matrix(data = NA, nrow = nrow(chosen_temp), ncol = 19)
  new_matrix[,1] <-chosen_temp$year
  new_matrix[,2] <- chosen_temp$julian_day
  new_matrix[,3] <- chosen_temp$month
  new_matrix[,4] <- chosen_temp$SST
  new_matrix[,5] <- chosen_temp$tdiff
  new_matrix[,6] <- unlist(residence_depth2(SST = new_matrix[,4], new_matrix[,5]))
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
                                   WF = new_matrix[i,9]) # auf aktuelles Gewicht bezogen
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



###### example 1: ----
#### run example 1

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



###### example 2: numeric growth for one cod age 3 - 2016 to age 4 in 2017 ----

#read in length at age data set
smalk_data <- read.table("length_weight_at_age_females_13_19.csv",header = TRUE, sep = ";")

chosen_year <- 2016

end_length_vector <- rnorm(n = 25, mean = smalk_data$mean_length[smalk_data$Year == chosen_year+1 & smalk_data$Age == 4],
                      sd = smalk_data$sd_length[smalk_data$Year == chosen_year+1 & smalk_data$Age == 4])
start_length <- smalk_data$mean_length[smalk_data$Year == chosen_year & smalk_data$Age == 3]


start_weight <- ps_lwr_a *start_length ^ ps_lwr_b
phys_weight <- ps_lwr_a *start_length ^ ps_lwr_b

chosen_temp <- temp_raw_dat(y = chosen_year, raw_temp_dat)

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
new_matrix[1,8] <- start_length # start length
new_matrix[1,9] <- start_weight # start weight
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

plot(1:365,new_matrix[,11], type = "l") # stomach content [g]
lines(1:365,new_matrix[,13], col = 2) # daily consumption [g]
lines(1:365,new_matrix[,14], col =3) # daily maintenance ration [g]

#temp und scw
plot(1:365,new_matrix[,13], type = "l")
lines(1:365, new_matrix[,7], col = 2)
lines(1:365, new_matrix[,11], col =4)


###example 3 ----
smalk_data <- read.table("smalk_raw_data_females.csv",header = TRUE, sep = ";")

chosen_year <- 2013
chosen_age <- 3

end_length_vector <- sample(smalk_data$length[smalk_data$year == chosen_year+1 & smalk_data$age == chosen_age+1]/10, 100, replace = TRUE)
hist(end_length_vector)

start_length_vector <- smalk_data$length[smalk_data$year == chosen_year & smalk_data$age == chosen_age]/10
hist(start_length_vector)

new_list <- vector(mode = "list", length = 100)
for(i in 1:100){
  start_length <- sample(start_length_vector,1,replace = TRUE) # random start length
  test1 <- grow(start_length = start_length, chosen_year = chosen_year)
  new_list[[i]] <-test1
}


#bind all matrices

j <- do.call("rbind",new_list)

j_data <- as.data.frame(j)
names(j_data) <- c("year","julian_day","month","SST","tdiff","res_depth","temp_at_depth","start_length",
                   "start_weight","phys_weight","scw","diet_clust","consumption","rmaint","ffg","weight_gain","new_weight",
                   "pred_length","growth_day")
summary(j_data)
library(plyr)
library(ggplot2)

#check <- ddply(j_data, .(growth_day), summarise, mean_resdepth= mean(res_depth), sd_resdepth = sd(res_depth))

#residence_depth_plot <-ggplot(data = check, aes(x = growth_day, y = mean_resdepth))+geom_line()+
#  geom_ribbon(aes(ymin=check$mean_resdepth - check$sd_resdepth, 
#                  ymax=check$mean_resdepth + check$sd_resdepth), alpha=.2)+theme_bw()+
#  ylab("Residence depth [m]") + xlab("modelling time [days]") + ylim(35,0)

###calculate it with 5 and 95% percentile instead of standard deviation
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


cowplot::plot_grid(residence_depth_plot2, rel_cons_plot , rel_maint_plot, 
                   rel_ffg_plot, boxplot_end_length, labels = "auto")

shapiro.test(end_length_vector) 
shapiro.test(j_data$pred_length[j_data$julian_day == 45])

wilcox.test(j_data$pred_length[j_data$julian_day == 45], end_length_vector)

plot(j_data$growth_day,j_data$new_weight)

#### grow and plotting as function ----
smalk_data <- read.table("smalk_raw_data_females.csv",header = TRUE, sep = ";")

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
  library("plyr")
  library("ggplot2")
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
  
  ##add some statistical tests
  #test on gaussian distribution of pred. and obs. lengths
  shap_pred_length <- shapiro.test(floor(j_data$pred_length[j_data$julian_day == 45]))[2]
  shap_obs_length <-shapiro.test(end_length_vector)[2]
  #test on homogeneity of variances across groups
  dat <- rbind(data.frame(length = floor(j_data$pred_length[j_data$julian_day == 45]), 
                         group = "predicted"),
              data.frame(length = end_length_vector,  group = "observed"))
  library("car")
  levene_obs_pred <- leveneTest(y =dat$length, group = dat$group)[3]$`Pr(>F)`[1]
  t_test_obs_pred <- t.test(x = floor(j_data$pred_length[j_data$julian_day == 45]), 
                            y = end_length_vector)[3]
  mwu_test_obs_pred <- wilcox.test(x = floor(j_data$pred_length[j_data$julian_day == 45]), 
                                   y = end_length_vector)[3]
  new_tibble <- tibble::tibble(list(j_data), pred_plot, list(end_length_vector),
         shap_pred_length, shap_obs_length,
         levene_obs_pred, t_test_obs_pred,
        mwu_test_obs_pred, n_start_indiv, n_end_indiv,
        chosen_age, chosen_year, number_ind)
  return(new_tibble)
  
}
test <- grow_and_plot(chosen_year = 2016, chosen_age = 3, number_ind = 25,
                      smalk_data = smalk_data)


###example 4: run models for males ----
smalk_data <- read.table("smalk_raw_data_males.csv",header = TRUE, sep = ";")

## use these coefficients to calculate start weight (After spawning!)
ps_lwr_a <- 0.01052065
ps_lwr_b <- 2.96939058

## length weight relationship coefficients calculated with smalk data (only males) for spawners (s)
## over Sd22-23, 1991-2019 (R^2 = 0.950)
s_lwr_a <- 0.007808755
s_lwr_b <- 3.049125885


test <- grow_and_plot(chosen_year = 2015, chosen_age = 3, number_ind = 10,
                      smalk_data = smalk_data)
saveRDS(test, file = "year2016_age3_ind100c_males.rds")


test <- readRDS(file = "year2016_age3_ind100b_males.rds")
test$pred_plot # prediction_plot_grid
test$shap_pred_length[1] # shapiro test p-value  for predicted end lengths
test$shap_obs_length[1] # shapiro test p-value for observed end lengths
test$levene_obs_pred[1] # levene test for test on variance homogeneity among groups
test$t_test_obs_pred[1] # t test p-value
test$mwu_test_obs_pred[1] # mwu test p-value
test$n_start_indiv[1] # number of individuals observed in start length
test$n_end_indiv[1] # number of individuals observed in end length
test$number_ind[1] # number of individuals used in grwoth modelling
cowplot::ggsave(plot = test$pred_plot, filename = "age4_2014_males.jpg")


### run model for males in years 2013 to 2018
### age 2 to age 5
### for 100 individuals

for(i in 2013: 2017){
  for(k in 2:4){
    #male lwr parameters
    ps_lwr_a <- 0.01052065
    ps_lwr_b <- 2.96939058
    s_lwr_a <- 0.007808755
    s_lwr_b <- 3.049125885
    # numer of modelled individuals
    n_ind <- 100
    # start modelling
    test <- grow_and_plot(chosen_year = i, chosen_age = k, number_ind = n_ind,
                          smalk_data = smalk_data)
    saveRDS(test, file = paste("age",k,"year",i,"ind",n_ind,".rds" ,sep = "_"))
    ggsave(plot = test$pred_plot, filename = paste("age",k,"year",i,"ind",n_ind,".png" ,sep = "_"))
  }
}

