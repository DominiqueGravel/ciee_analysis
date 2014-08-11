
original_data = read.csv2("data/marinesize.csv",dec=".")
MPred = log10(original_data$predator_length)
MPrey = log10(original_data$prey_length)

# Filter out strange prey lengths
test = numeric(length(MPrey))
test[MPrey>1.8 & MPred < 1.5] = 1
MPred = MPred[test==0]
MPrey = MPrey[test==0]


#### 
# Without temperature effects
Pars = read.table("data/ParsLengthNoTemp.txt")
seqM = seq(min(MPrey),max(MPred),0.01)
l1 = Pars$a0 + Pars$a1*seqM
l2 = Pars$a0 + Pars$a1*seqM + 3*(Pars$b0 + Pars$b1*seqM)
l3 = Pars$a0 + Pars$a1*seqM - 3*(Pars$b0 - Pars$b1*seqM)

quartz(height = 6, width = 7)
par(mar=c(6,6,2,1))	
plot(MPred,MPrey,pch=21,cex.axis = 1.25, cex.lab = 1.5, xlab = "Predator size", ylab = "Prey size")
lines(seqM,l1, col = "darkred", lwd = 2)
lines(seqM,l2, col = "gray", lwd = 2)
lines(seqM,l3, col = "gray", lwd = 2)
dev.copy2pdf(file = "figures/PredPreyNoTemp.pdf")

#### 
# With temperature effects
Pars = read.table("data/ParsLengthTemp.txt")
seqM = seq(min(MPrey),max(MPred),0.01)
l1 = Pars$a0 + (Pars$a1 + Pars$a2*25)*seqM + Pars$a3*25
l2 = Pars$a0 + (Pars$a1 + Pars$a2*10)*seqM + Pars$a3*10

#Temp = original_data$mean_annual_temp
#RK = rank(Temp)
vec.col = rainbow(n = length(Temp), alpha = 0.1, start = 0.3, end = 1)[RK]

quartz(height = 6, width = 7)
par(mar=c(6,6,2,1))	
plot(MPred,MPrey,pch=19,cex.axis = 1.25, cex.lab = 1.5, xlab = "Predator size", ylab = "Prey size",col = rgb(0, 0, 0, alpha = 0.25),cex = 0.2)
lines(seqM,l1, col = "darkred", lwd = 2)
lines(seqM,l2, col = "darkblue", lwd = 2)
legend("topleft",lty = 1, col = c("darkred","darkblue"), legend = c("T = 25", "T = 15"), lwd = 2, bty = "n")

dev.copy2pdf(file = "figures/PredPreyTemp.pdf")


