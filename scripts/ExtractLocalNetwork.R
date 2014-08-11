
# Load spatial data
load("data/Data_SST_1961_1980_biomod_plato")
env = Data_1961_1980_biomod_plato
Temp = (env$Summer*3 + env$October_Juin*9)/12

# Load size data
Size = read.table("data/sub_size_med.txt", header = T)

# Load presence-absence data
Pres = read.table("data/sub_pres_med.txt")

# Load parameters
Pars = read.table("data/ParsLengthTemp.txt")

###
ExtractLocalNetwork = function(X,Temp,Pars,Pres,Size) {

	# Extract the local species list
	LocalPres = Pres[X,]
	S = sum(LocalPres)
	
	# Extract body size for this site
	LocalM = Size[LocalPres == 1,3]

	# Compute the web for this particular temperature
	MPrey = log10(expand.grid(LocalM,LocalM)[,1])
	MPred = log10(expand.grid(LocalM,LocalM)[,2])

	# Compute interaction probability
	pL = pLFitted(MPrey,MPred,Pars,Temp = Temp[X])

	# Transform into an interaction matrix 
	L = matrix(pL, nr = S, nc = S, byrow = FALSE)

	return(L)
}

WEBS = list()
for(i in 1:8154) WEBS[[i]] = ExtractLocalNetwork(i,Temp,Pars,Pres,Size)








