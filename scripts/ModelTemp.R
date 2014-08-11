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


######################
# The model
Model = function(MPrey,MPred,meanMprey,sdMprey,a0,a1,a2,a3,b0,b1,h,Temp) {

	# Optimum and range
#	o = a0 + a1*MPred 
	o = a0 +(a1 + a2*Temp)*MPred + a3*Temp 
	r = b0 + b1*MPred
		
	# Compute the conditional
	pLM = h*exp(-(o-MPrey)^2/2/r^2)

	# Compute the marginal
	pM = dnorm(x=MPrey,mean=meanMprey,sd=sdMprey)
	
	#  Integrate the denominator
	pL = r/(r^2+sdMprey^2)^0.5*h*exp(-(o-meanMprey)^2/2/(r^2+sdMprey^2))	
	
	# Compute the posterior probability
	pML = pLM*pM/pL
	
	return(pML)		
}
######################

######################
PDF = function(pML) log(pML)
######################

######################
pLFitted = function(MPrey,MPred,Pars,Temp) {
	with(Pars, {
	# Optimum and range
	o = a0 +(a1 + a2*Temp)*MPred + a3*Temp 
#	o = a0 + a1*MPred 
	r = b0 + b1*MPred 
	
	# Compute the conditional
	h*exp(-(o-MPrey)^2/2/r^2)
	})
}
######################
