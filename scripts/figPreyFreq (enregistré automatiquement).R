
original_data = read.csv2("data/size_Barnes2008.csv",dec=".")
MPred = log10(original_data$predator_length)
MPrey = log10(original_data$prey_length)
quartz(height = 6, width = 7)
par(mar=c(5,6,2,1))	
hist(MPrey,xlab = expression(M["prey"]),main="",cex.axis = 1.25, cex.lab = 1.5,freq = FALSE)
M = seq(-1,4,0.01)
P = 0.4*exp(-(0-M)^2/2/0.5^2)
lines(M,P,col = "darkred",lwd = 2)
dev.copy2pdf(file = "figures/prey_freq.pdf")

M = c(unique(log10(original_data$predator_length)),unique(log10(original_data$prey_length)))


MPrey = expand.grid(M,M)[,1]
MPred = expand.grid(M,M)[,2]
PreyNames = expand.grid(Size[,1],Size[,1])[,1]
PredNames = expand.grid(Size[,1],Size[,1])[,2]




