load("data/Data_SST_1961_1980_biomod_plato")
env = Data_1961_1980_biomod_plato
Temp = (env$Summer*3 + env$October_Juin*9)/12

# Load presence-absence data
Pres = read.table("data/sub_pres_med.txt")
S = ncol(Pres)

library(randomForest)
Pres_pred = Pres*0

for(i in 1:S) {
	set.seed(23)
	data = data.frame(pres = Pres[,i],Env = env[,4:10])
	SDM = randomForest(pres ~ ., data, ntree = 100)
	set.seed(23)
	Pres_pred[,i] = predict(SDM)
	cat(i,'\n')
}

write.table(Pres_pred, file = "data/PresPredMed.txt")









