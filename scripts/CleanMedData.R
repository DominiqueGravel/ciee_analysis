source("scripts/model.R")
size_med = read.table("data/size_med.txt", header = T,row.names = 1)
size_med = cbind(row.names(size_med)[order(size_med)],size_med[order(size_med),1])
sp_list = read.csv2("data/sp_list.csv")[,1:2]

# Load presence-absence data
pres = read.table("data/pres_med.txt")

# Match the column names with body size
ESP = as.character(sp_list[match(size_med[,1],as.character(sp_list[,1])),2])

# Remove species with no match
sub_size_med = subset(size_med,ESP!="NA") 
sub_ESP = subset(ESP,ESP!="NA")
sub_pres = pres[,match(sub_ESP,names(pres))]

# Record clean files
write.table(cbind(sub_ESP,sub_size_med), "data/sub_size_med.txt")
write.table(sub_pres, "data/sub_pres_med.txt")


