# Load size data
Size = read.table("data/sub_size_med.txt", header = T)

# Load presence-absence data
Pres = read.table("data/pres_med.txt")

source("scripts/Model.R")

# Expand all pairs of species 
M = as.numeric(Size[,3])
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
	sum(Pres[,col1]*Pres[,col2]) 
}
Cooc = apply(cbind(as.character(PreyNames),as.character(PredNames)),1,wrapper)
pCooc = Cooc/8154

pCoocSDM = read.table("data/pCoocSDM.txt")/8154

# Report the results
Niche = read.table("data/SpMeanTemp.txt")
PreyNiche=expand.grid(Niche[,1],Niche[,1])[,1]
PredNiche=expand.grid(Niche[,1],Niche[,1])[,2]
res = cbind(as.character(PreyNames),as.character(PredNames),pL,pCooc,pLCoocSDM,PreyNiche,PredNiche)

write.table(res,"data/metaweb.txt")