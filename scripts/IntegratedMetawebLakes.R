# Load size data
Size = read.table("data/SizeLakes.txt", header = T)

# Load presence-absence data
Pres = read.table("data/PresLakes.txt",header = T)

source("scripts/ModelNoTemp.R")

# Expand all pairs of species 
M = log10(as.numeric(Size[,3]))
MPrey = expand.grid(M,M)[,1]
MPred = expand.grid(M,M)[,2]
PreyNames = expand.grid(Size[,1],Size[,1])[,1]
PredNames = expand.grid(Size[,1],Size[,1])[,2]

# Compute interaction probability
Pars = read.table("data/ParsLengthNoTemp.txt")
Pars$a2 = 0
Pars$a3 = 0
pL = pLFitted(MPrey,MPred,Pars,Temp = 30)

# Compute co-occurrence
wrapper = function(x){
	col1 = match(x[1],names(Pres))
	col2 = match(x[2],names(Pres))
	if(is.na(col1)==FALSE & is.na(col2)==FALSE) sum(Pres[,col1]*Pres[,col2]) 
	else NA
}
Cooc = apply(cbind(as.character(PreyNames),as.character(PredNames)),1,wrapper)
pCooc = Cooc/nrow(Pres)

# Report the results

#Niche = read.table("data/SpMeanTemp.txt")
#PreyNiche=expand.grid(Niche[,1],Niche[,1])[,1]
#PredNiche=expand.grid(Niche[,1],Niche[,1])[,2]
res = cbind(as.character(PreyNames),as.character(PredNames),pL,pCooc)

write.table(res,"data/metawebLakes.txt")

S = length(pL)^0.5
L = matrix(pL,nr = S, nc = S, byrow = FALSE)
image(c(1:S),c(1:S),L)




