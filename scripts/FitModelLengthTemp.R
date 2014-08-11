######################
# Source packages
rm(list = ls())
source("scripts/anneal.R")
source("scripts/likeli.R")
source("scripts/analyze_function.R")
source("scripts/likdisplay.R")
source("scripts/predicted_results.R")
source('scripts/Model.R')

# Load data
original_data = read.csv2("data/size_Barnes2008.csv",dec=".")
MPred = log10(original_data$predator_length)
MPrey = log10(original_data$prey_length)
Temp = original_data$mean_annual_temp

# Filter out strange prey lengths
test = numeric(length(MPrey))
test[MPrey>1.8 & MPred < 1.5] = 1
MPred = MPred[test==0]
MPrey = MPrey[test==0]
Temp = Temp[test==0]

data = data.frame(MPrey=MPrey,MPred=MPred, Temp=Temp)

# Initial values 
lm_M = lm(MPrey~MPred)
par = list()
par$a0 = lm_M$coefficients[1]
par$a1 = lm_M$coefficients[2]
par$a2 = 0
par$a3 = 0
par$b0 = sd(lm_M$residuals)
par$b1 = 0
par$h = 1

# Boundaries
par_lo = list(a0 = -100, a1 = 0, a2 = -10, a3 = -10, b0 = 0, b1 = -10, h = 0)
par_hi = list(a0 = 100, a1 = 10, a2 = 10, a3 = 10, b0 = 100, b1 = 10, h = 1)

#par_lo = list(a0 = -100, a1 = 0, b0 = 0, b1 = -10, h = 0)
#par_hi = list(a0 = 100, a1 = 10, b0 = 100, b1 = 10, h = 1)

# Stuff for the anneal function
var = list()
var$MPrey = "MPrey"
var$MPred = "MPred" 
var$pML = "predicted" 
var$Temp = "Temp"
var$meanMprey = mean(MPrey)
var$sdMprey = sd(MPrey)

# Maximum likelihood estimation
test = anneal(model = Model, par = par, var = var, source_data = data, 
par_lo = par_lo, par_hi = par_hi, dep_var = "MPrey", pdf = PDF, 
max_iter = 10000, hessian = FALSE, initial_temp = 0.1)
write.table(test$best_pars,file = "data/ParsLengthTemp.txt")



