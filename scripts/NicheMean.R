
# Load spatial data
load("data/Data_SST_1961_1980_biomod_plato")
env = Data_1961_1980_biomod_plato
Temp = (env$Summer*3 + env$October_Juin*9)/12

# Load size data
Size = read.table("data/sub_size_med.txt", header = T)

# Load presence-absence data
Pres = read.table("data/sub_pres_med.txt")

SpMeanTemp = apply(Temp*Pres,2,sum)/apply(Pres,2,sum)

write.table(SpMeanTemp,"data/SpMeanTemp.txt")



