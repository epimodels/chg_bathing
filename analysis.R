###################################################################
# Statistical Analysis of CHG Efficacy in ICU MRSA Decolonization #
###################################################################

## Import Packages and Import Data
library(vioplot)

## Main Estimates and Plotting for CHG and Mupirocin Efficacy

# Delta and Zeta estimates for main and latent model
#delta <- as.data.frame(read_csv("delta_main.csv",header=FALSE,col_names=c("Estimate")))

delta <- read.table("deltamain.csv",header=FALSE,col.names=c("Estimate"))
delta_l <-read.table("deltamain_latent.csv",header=FALSE,col.names=c("Estimate"))
zeta <- read.table("zetamain.csv",header=FALSE,col.names=c("Estimate"))
zeta_l <- read.table("zetamain_latent.csv",header=FALSE,col.names=c("Estimate"))

# Frequency sensitivity analysis
frequency <- as.data.frame(read.csv(("frequency.csv")))
frequency$FactorFreq <- as.character(frequency$Frequency)

# Parameter Estimates for Delta and Zeta
delta_est <- quantile(delta$Estimate, prob = c(0.025, 0.50, 0.975))
delta_latent_est <- quantile(delta_l$Estimate, prob = c(0.025, 0.50, 0.975))
zeta_est <- quantile(zeta$Estimate, prob = c(0.025, 0.50, 0.975))
zeta_latent_est <- quantile(zeta_l$Estimate, prob = c(0.025, 0.50, 0.975))

# Plots of Posteriors

par(mfrow=c(2,2))

delta_den <- density(delta$Estimate)
delta_l_den <- density(delta_l$Estimate)
delta_y_max <- max(delta_den$y,delta_l_den$y)

zeta_den <- density(zeta$Estimate)
zeta_l_den <- density(zeta_l$Estimate)
zeta_y_max <- max(zeta_den$y,zeta_l_den$y)

plot(delta_den,lwd=3,main="CHG Posterior - Main Model",xlab="Per-Use CHG Efficacy",xlim=c(0,1),cex.lab=1.5,ylim=c(0,delta_y_max))
abline(v=median(delta$Estimate), col=c("red3"),lwd=3,lty=3)
legend("topright", c("Posterior Density","Median"), lwd=3, col=c("black","red3"),lty=c(1,3), bty='n', cex=1)

plot(delta_l_den,lwd=3,main="CHG Posterior - Latency Model",xlab="Per-Use CHG Efficacy",xlim=c(0,1),cex.lab=1.5,ylim=c(0,delta_y_max))
abline(v=median(delta_l$Estimate), col=c("red3"),lwd=3,lty=3)
legend("topright", c("Posterior Density","Median"), lwd=3, col=c("black","red3"),lty=c(1,3), bty='n', cex=1)

plot(zeta_den,lwd=3,main="Mupirocin Posterior - Main Model",xlab="Per-Use Mupirocin Efficacy",xlim=c(0,1),cex.lab=1.5,ylim=c(0,zeta_y_max))
abline(v=median(zeta$Estimate), col=c("red3"),lwd=3,lty=3)
legend("topright", c("Posterior Density","Median"), lwd=3, col=c("black","red3"),lty=c(1,3), bty='n', cex=1)

plot(zeta_l_den,lwd=3,main="Mupirocin - Latency Model",xlab="Per-Use Per-Use Mupirocin Efficacy",xlim=c(0,1),cex.lab=1.5,ylim=c(0,zeta_y_max))
abline(v=median(zeta_l$Estimate), col=c("red3"),lwd=3,lty=3)
legend("topright", c("Posterior Density","Median"), lwd=3, col=c("black","red3"),lty=c(1,3), bty='n', cex=1)

## Exploration of Sensitivity of Results to Timing of Application
# ANOVA and Tukey HSD for the differences in the main scenarios
AcquisitionTest <- aov(frequency$Cases ~ frequency$FactorFreq)

print(summary(AcquisitionTest))
print(TukeyHSD(AcquisitionTest))

control <- subset(frequency,FactorFreq=="0")
daily <- subset(frequency, FactorFreq=="24")
bidaily <- subset(frequency, FactorFreq=="48")
threeday <- subset(frequency,FactorFreq=="72")
fourday <- subset(frequency,FactorFreq=="96")
fiveday <- subset(frequency,FactorFreq=="120")

par(mfrow=c(1,1))
vioplot(control$Cases,daily$Cases,bidaily$Cases,threeday$Cases,fourday$Cases,fiveday$Cases,names=c("Baseline","24","48","72","96","120"),
        col="Grey90",drawRect=FALSE,ylim=c(0,3))
title(ylab="MRSA Acquisitions (per 1,000 patient-days)",cex.lab=1.5)
title(xlab="Time Between Bathing (Hours)",cex.lab=1.5)
segments(0.75,mean(control$Cases),1.25,mean(control$Cases),col="black",lwd=3)
segments(1.75,mean(daily$Cases),2.25,mean(daily$Cases),col="black",lwd=3)
segments(2.75,mean(bidaily$Cases),3.25,mean(bidaily$Cases),col="black",lwd=3)
segments(3.75,mean(threeday$Cases),4.25,mean(threeday$Cases),col="black",lwd=3)
segments(4.75,mean(fourday$Cases),5.25,mean(fourday$Cases),col="black",lwd=3)
segments(5.75,mean(fiveday$Cases),6.25,mean(fiveday$Cases),col="black",lwd=3)

legend("topright", c("Scenario Mean"), lwd=3,col=c("black"),lty=1, bty='n', cex=1.25)

## Global Sensitivitity Analysis
# Global parameter sensitivity analysis
global_headers = c("Estimate","Rho_N","Rho_D","Psi","Theta","Nu","Iota_N","Iota_D","Tau_N","Tau_D","Mu")
global <- read.csv("sensweep.csv",header=FALSE,col.names=global_headers)
global <- global[-which(global$Estimate == 0 ), ] #remove zeroes - runs where no estimate change was found

# Normalize Global Sensitivity to $-change
global$Rho_N<- (global$Rho_N-1)*100
global$Rho_D<- (global$Rho_D-1)*100
global$Psi<- (global$Psi-1)*100
global$Theta<- (global$Theta-1)*100
global$Nu<- (global$Nu-1)*100
global$Iota_N<- (global$Iota_N-1)*100
global$Iota_D<- (global$Iota_D-1)*100
global$Tau_D<- (global$Tau_D-1)*100
global$Tau_N<- (global$Tau_N-1)*100
global$Mu<- (global$Mu-1)*100

par(mfrow=c(1,1))
sens_model <- lm(Estimate ~ Rho_N + Rho_D + Psi + Theta + Nu + Iota_N + Iota_D + Tau_N + Tau_D + Mu,data=global)
summary(sens_model)

# Plot Sensitivity
coefficients <- sort(sens_model$coefficients[2:11]*10)
cols <- c("grey90", "grey50")[(coefficients < 0) + 1]  
greek_names1 <- c(expression(iota[N]),expression(tau[N]),expression(iota[D]),
                  expression(mu),expression(tau[D]),expression(rho[D]),
                  expression(theta),expression(nu),expression(psi),expression(rho[N]))

barplot(coefficients, horiz=TRUE,names.arg=greek_names1,ylab="Parameter",xlab="Change in Combined Estimate",
        cex.lab=1.25,cex.axis=1.25,cex.names=1.75,col=cols,xlim=c(-0.05,0.10))