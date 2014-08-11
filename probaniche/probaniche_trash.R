######################
# Generate fake data
rm(list = ls())
M = expand.grid(seq(0,100,1),seq(0,100,1))
pred = as.matrix(M[,1])
prey = as.matrix(M[,2])
h = 0.5
a0 = -10
a1 = 0.75
b0 = 10
b1 = 0.0

o = a0 + a1*pred
r = b0 + b1*pred
pLM = h*exp(-(o-prey)^2/2/r^2)
L = numeric(nrow(M))
rand = runif(nrow(M),0,1)
L[pLM > rand] = 1

M = M[L==1,]
meanMprey = mean(M[,1])
sdMprey = sd(M[,1])
par(mar=c(5,5,2,1))
plot(M[,2],M[,1],pch=19,ylim = c(0,100),xlim = c(0,100))
abline(0,1)

######################
# The conditional probability function
pLM_fn = function(Mprey,Mpred,a0,a1,b0,b1,h) {
	o = a0 + a1*Mpred
#	r = b0 + b1*Mpred
	r = b0 
	h*dnorm(x = Mprey, mean = o, sd = r)
}

######################
# Denominator
denom = function(x,Mpred,meanMprey,sdMprey,a0,a1,b0,b1,h) pLM_fn(Mprey=x,Mpred=Mpred,a0=a0,a1=a1,b0=b0,b1=b1,h=h)*dnorm(x=x,mean=meanMprey,sd=sdMprey)

######################
# The marginal probability function
#pL_fn = function(Mpred,meanMprey,sdMprey,a0,a1,b0,b1,h) integrate(f=denom,lower=-Inf,upper=Inf,Mpred=Mpred,a0=a0,a1=a1,b0=b0,b1=b1,h=h,meanMprey=meanMprey,sdMprey=sdMprey)[1]	

######################
# The model
# model = function(X,prey,pred,meanMprey,sdMprey,a0,a1,b0,b1,h) {

	# Mprey = prey[X]
	# Mpred = pred[X]

	# # Compute the conditional
	# pLM = pLM_fn(Mprey=Mprey,Mpred=Mpred,a0=a0,a1=a1,b0=b0,b1=b1,h=h)

	# # Compute the marginal
	# pM = dnorm(x=Mprey, mean = meanMprey,sd = sdMprey)
	
	# #  Integrate the denominator
	# pL = as.numeric(pL_fn(Mpred=Mpred,meanMprey=meanMprey,sdMprey=sdMprey,a0=a0,a1=a1,b0=b0,b1=b1,h=h))
	
	# # Compute the posterior probability
	# pML = pLM*pM/pL
	
	# return(pML)		
# }

# wrapper = function(prey,pred,meanMprey,sdMprey,a0,a1,b0,b1,h) {
	# X = as.matrix(c(1:length(prey)),nc = 1)
	# apply(X,MARGIN=1,FUN=model,prey=prey, pred=pred,meanMprey=meanMprey,sdMprey=sdMprey,a0=0,a1=a1,b0=b0,b1=b1,h=h)
# }

model = function(Mprey,Mpred,meanMprey,sdMprey,a0,a1,b0,b1,h) {

	# Optimum and range
	o = a0 + a1*Mpred
	r = b0 	
	
	# Compute the conditional
	pLM = h*dnorm(x = Mprey, mean = o, sd = r)
	
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
alldata = read.table("data.txt",header=T)

# Initial values and controls
par = list()

# Compute the pred-prey relationship

# Compute the sd of the residuals

par$a0 = a0
par$a1 = a1
par$b0 = b0
par$b1 = b1
par$h = h
par_lo = list(a0 = -100, a1 = 0, b0 = 1, b1 = 0, h = 0)
par_hi = list(a0 = 100, a1 = 10, b0 = 100, b1 = 10, h = 1)

var = list()
var$Mprey = "Mprey"
var$Mpred = "Mpred" 
var$pML = "predicted" 
var$meanMprey = mean(M[,1])
var$sdMprey = sd(M[,1])

Mprey = M[,1]
Mpred = M[,2]
data = data.frame(Mprey = Mprey, Mpred = Mpred)

# Maximum likelihood estimation
test = anneal(model = model, par = par, var = var, source_data = data, 
par_lo = par_lo, par_hi = par_hi, dep_var = "Mprey", pdf = PDF, 
max_iter = 10000, hessian = FALSE, initial_temp = 1)
	


