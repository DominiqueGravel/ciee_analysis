
source("scripts/Model.R")

# Load size data
Size = read.table("data/sub_size_med.txt", header = T)

# Expand all pairs of species 
M = as.numeric(Size[,3])
MPrey = expand.grid(M,M)[,1]
MPred = expand.grid(M,M)[,2]
S = length(M)

Pars = read.table("data/ParsLengthTemp.txt")

Tp = seq(5,30,0.1)
res = numeric(length(Tp))
for(i in 1:length(Tp)) res[i] = sum(pLFitted(MPrey,MPred,Pars,Temp = Tp[i]))/S^2

quartz(height = 6, width = 7)
par(mar=c(5,6,2,1))	
plot(Tp,res,xlab = "Temperature", ylab = "Connectance",type = "l", cex.axis = 1.25, cex.lab = 1.5,lwd = 2)
dev.copy2pdf(file = "figures/ConnectanceTemp.pdf")
	

