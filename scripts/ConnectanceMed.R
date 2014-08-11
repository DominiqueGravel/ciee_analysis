# Environmental data
load("data/Data_SST_1961_1980_biomod_plato")
env = Data_1961_1980_biomod_plato
Temp = (env$Summer*3 + env$October_Juin*9)/12

######################
# Expected distribution data
Pres_pred = read.table("data/PresPredMed.txt")
expS = apply(Pres_pred,1,sum)
S = ncol(Pres_pred)

######################
# Expected interactions
source("scripts/Model.R")
Size = read.table("data/sub_size_med.txt", header = T)
Pars = read.table("data/ParsLengthTemp.txt")
M = as.numeric(Size[,3])
MPrey = expand.grid(M,M)[,1]
MPred = expand.grid(M,M)[,2]
PreyNames = expand.grid(Size[,1],Size[,1])[,1]
PredNames = expand.grid(Size[,1],Size[,1])[,2]
S = length(M)

# Save the metawebs at low and high temp
pLLow = pLFitted(MPrey,MPred,Pars,Temp = 10)
pLHigh = pLFitted(MPrey,MPred,Pars,Temp = 25)
mwLow = subset(cbind(PreyNames,PredNames,pLLow),pLLow>0.4)
mwHigh = subset(cbind(PreyNames,PredNames,pLHigh),pLHigh>0.4)
nrow(mwLow)/65536
nrow(mwHigh)/65536
write.table(mwLow, file = "data/MetawebLowTemp.txt")
write.table(mwHigh, file = "data/MetawebHighTemp.txt")

# Connectance without considering the effect of distribution
Tp = seq(min(Temp),max(Temp),0.01)
C = numeric(length(Tp))
for(i in 1:length(Tp)) C[i] = sum(pLFitted(MPrey,MPred,Pars,Temp = Tp[i]))/S^2

# Account for co-occurrence in each location
totalCooc = matrix(0,nr = S, nc = S)
res = matrix(nr = 8154, nc = 5)
for(z in 1:8154) {	
	pL = pLFitted(MPrey,MPred,Pars,Temp = Temp[z])
	matPrey = matrix(as.numeric(Pres_pred[z,]),nr = S, nc = S,byrow = FALSE)
	matPred = matrix(as.numeric(Pres_pred[z,]),nr = S, nc = S,byrow = TRUE)
	totalCooc = totalCooc + matPrey*matPred
	pCooc = stack(as.data.frame(matPrey*matPred))[,1]
	pObs = pL*pCooc
	res[z,1] = sum(Pres_pred[z,])
	res[z,2] = sum(pL)
	res[z,3] = sum(pCooc)
	res[z,4] = sum(pObs)
	res[z,5] = sum(pObs)/sum(Pres_pred[z,])^2
	cat(z,'\n')
}
write.table(totalCooc, file = "data/pCoocSDM.txt")
write.table(res,file = "data/ConnectanceMed.txt")


######################
# Plot the results
quartz(height = 6, width = 7)
par(mar=c(5,6,2,1))	
plot(Tp,C,xlab = "Temperature", ylab = "Connectance",type = "l", cex.axis = 1.25, cex.lab = 1.5,lwd = 2,ylim = c(0.14,0.29))
points(Temp,res[,5],pch = 21,bg = "gray")
dev.copy2pdf(file = "figures/ConnectanceTempSDM.pdf")

quartz(height = 6, width = 7)
par(mar=c(5,6,2,1))	
plot(Tp,C,xlab = "Temperature", ylab = "Connectance",type = "l", cex.axis = 1.25, cex.lab = 1.5,lwd = 2,ylim = c(0.23,0.29))
dev.copy2pdf(file = "figures/ConnectanceTemp.pdf")
	


# Return the global metaweb

