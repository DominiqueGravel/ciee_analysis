############################################
#
# Bayesien model of trophic interactions
# based on body size, using observation of interactions 
# only (presence)
# 
# Dominique Gravel
# June 25th, 2014
#
############################################


# source('model.R')
# source('web_construction.R')


######################
# The model
model = function(Mprey,Mpred,meanMprey,sdMprey,a0,a1,a2,a3,b0,b1,h,Temp) {

	# Optimum and range
	o = a0 + a1*Mpred 
	o = a0 +(a1 + a2*Temp)*Mpred + a3*Temp 
#	r = b0
	r = b0 + b1*Mpred
		
	# Compute the conditional
	pLM = h*exp(-(o-Mprey)^2/2/r^2)

	# Compute the marginal
	pM = dnorm(x=Mprey,mean=meanMprey,sd=sdMprey)
	
	#  Integrate the denominator
	pL = r/(r^2+sdMprey^2)^0.5*h*exp(-(o-meanMprey)^2/2/(r^2+sdMprey^2))	
	
	# Compute the posterior probability
	pML = pLM*pM/pL
	
	return(pML)		
}

# PDF
PDF = function(pML) log(pML)

######################
# Source packages
setwd("/Users/DGravel/Desktop/probaniche")
source("anneal.R")
source("likeli.R")
source("analyze_function.R")
source("likdisplay.R")
source("predicted_results.R")

# Load data
original_data = read.csv2("marinesize.csv",dec=".")
Mpred = log10(original_data$predator_mass)
Mprey = log10(original_data$prey_mass)
Temp = original_data$mean_annual_temp
data = data.frame(Mprey=Mprey,Mpred=Mpred, Temp=Temp)

# Initial values 
lm_M = lm(Mprey~Mpred)
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

# Stuff for the anneal function
var = list()
var$Mprey = "Mprey"
var$Mpred = "Mpred" 
var$pML = "predicted" 
var$Temp = "Temp"
var$meanMprey = mean(Mprey)
var$sdMprey = sd(Mprey)

# Maximum likelihood estimation
test = anneal(model = model, par = par, var = var, source_data = data, 
par_lo = par_lo, par_hi = par_hi, dep_var = "Mprey", pdf = PDF, 
max_iter = 10000, hessian = FALSE, initial_temp = 1)
write.table(test$best_pars,file = "best_pars.txt")


######################
# Plot the results
quartz(height = 6, width = 7)
par(mar=c(6,6,2,1))	
plot(Mpred,Mprey,pch=21,cex.axis = 1.25, cex.lab = 1.5, xlab = "Predator size", ylab = "Prey size")
seqM = seq(min(Mpred),max(Mpred),0.001)
lines(seqM, test$best_pars$a0 + test$best_pars$a1*seqM, col = "darkred",lwd = 2)
lines(seqM, test$best_pars$a0 + test$best_pars$a1*seqM + 3*(test$best_pars$b0+test$best_pars$b1*seqM), col = "gray")
lines(seqM, test$best_pars$a0 + test$best_pars$a1*seqM - 3*(test$best_pars$b0+test$best_pars$b1*seqM), col = "gray")

dev.copy2eps(file = "PredPreyFit.pdf")

######################
get_opt = function(Mpred, best_pars,Temp) {
	a0 = best_pars$a0
	a1 = best_pars$a1
	a2 = best_pars$a2
	a3 = best_pars$a3		
	a0 +(a1 + a2*Temp)*Mpred + a3*Temp 
}

quartz(height = 6, width = 7)
par(mar=c(6,6,2,1))	
plot(Mpred,Mprey,pch=21,cex.axis = 1.25, cex.lab = 1.5, xlab = "Predator size", ylab = "Prey size")
lines(seqM,get_opt(seqM,test$best_pars,Temp= 15), col = "darkblue", lwd = 2)
lines(seqM,get_opt(seqM,test$best_pars,Temp= 30), col = "darkred", lwd = 2)
legend("topleft",lty = 1, col = c("darkred","darkblue"), legend = c("T = 30", "T = 15"), lwd = 2, bty = "n")
dev.copy2eps(file = "PredPreyTemp.pdf")

######################
# Reconstruct the probabilistic network

pL_fitted = function(Mprey,Mpred,a0,a1,a2,a3,b0,b1,h, Temp) {
	# Optimum and range
	o = a0 +(a1 + a2*Temp)*Mpred + a3*Temp 
	r = b0 + b1*Mpred 
	
	# Compute the conditional
	pLM = h*exp(-(o-Mprey)^2/2/r^2)
}

# Pair all potential species
a0 = test$best_pars$a0
a1 = test$best_pars$a1
a2 = test$best_pars$a2
a3 = test$best_pars$a3
b0 = test$best_pars$b0
b1 = test$best_pars$b1
h = test$best_pars$h

IDprey = log10(tapply(original_data$prey_mass,original_data$prey,mean))
IDpred = log10(tapply(original_data$predator_mass,original_data$predator,mean))
M = c(IDprey,IDpred)
sp_names = names(M)
S = length(M)

Mprey= rep(M, time = S)
Mpred = rep(M, each = S)

pL25 = pL_fitted(Mprey,Mpred,a0,a1,a2,a3,b0,b1,h,Temp = 25)
pL10 = pL_fitted(Mprey,Mpred,a0,a1,a2,a3,b0,b1,h,Temp = 10)

# Put a cutoff value

L25 = numeric(length(pL25))
L25[pL25 > 0.6] = 1
sum(L25)/S^2

L10 = numeric(length(pL10))
L10[pL10 > 0.6] = 1
sum(L10)/S^2


pairs25 = cbind(names(Mprey),names(Mpred))[L25==1,]
pairs10 = cbind(names(Mprey),names(Mpred))[L10==1,]

write.table(pairs25, file = "pairs25.txt")
write.table(pairs10, file = "pairs10.txt")




