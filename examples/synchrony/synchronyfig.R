## Started 16 August 2025 ##
## Taken from git/projects/misc/miscmisc/bayesianflowsexample/example.Rmd ## 

## housekeeping
rm(list=ls()) 
options(stringsAsFactors = FALSE)
set.seed(3777)

## libraries
library(rstan)
options(mc.cores = parallel::detectCores())

# setwd
setwd("/Users/lizzie/Documents/git/projects/misc/miscmisc/bayesianflows/examples/synchrony")

## flags
# You must set to true and RUN the models once for setting this to FALSE to work
runmodels <- FALSE

#######################
## Step 3: real data ##
#######################

# get the data
rawlong.tot2 <- read.csv("output/rawlong.tot2.csv")

# Formatting for R stan (several ways to do this, this is one)
N <- nrow(rawlong.tot2)
y <- rawlong.tot2$phenovalue
Nspp <- length(unique(rawlong.tot2$species)) #newid is character !
species <- as.numeric(as.factor(rawlong.tot2$species))
year <- rawlong.tot2$yr1981

if(runmodels){
# See the stan code on this model for notes on what it does
syncmodelhis <- stan("stan/twolevelhierslopeint.stan", data=c("N","Nspp","y","species","year"),
                   iter=4000, warmup=3000, chains=4, cores=4)
save(syncmodelhis, file="output/syncmodelhis.Rdata")
}

if(!runmodels){
load("output/syncmodelhis.Rdata")
}

##########################################
## Posterior predictive checks (Step 4) ##
##########################################

Nreal <- nrow(rawlong.tot2)
yreal <- rawlong.tot2$phenovalue


# First, plot the real data used in the model
# pdf("graphs/realdata_formodel.pdf", height=4, width=6)
par(mar=c(3,3,1,1), mgp=c(1.5,.5,0), tck=-.01)
plot(range(year), range(yreal), type="n", xlab="Year",
     ylab="Day of year", bty="l", main="Raw real data")
for (j in 1:Nspp){
  lines(year[species==j], yreal[species==j])
}
hist(yreal, xlab="Day of year", main="Real data")
# dev.off()

# What does a similar plot look like using the model output?
syncmodelhispost <- extract(syncmodelhis) 
# hist(syncmodelhispost$mu_b, xlab="Change per year")

# extract means for now (other ways to extract the mean)
sigma_y <- mean(syncmodelhispost$sigma_y) 
sigma_a <- mean(syncmodelhispost$sigma_a) 
sigma_b <- mean(syncmodelhispost$sigma_b) 
mu_b <- mean(syncmodelhispost$mu_b) 
mu_a <- mean(syncmodelhispost$mu_a) 

a <- rnorm(Nspp, mean=mu_a, sd=sigma_a)
b <- rnorm(Nspp, mean=mu_b, sd=sigma_b)

N <- Nreal

ypred <- length(N) 
for (n in 1:N){
    s <- species[n]
    ypred[n] <- a[s] + b[s]*year[n]
}
y <- rnorm(N, ypred, sigma_y)

pdf("graphs/rawvsonepredictivecheck.pdf", height=8, width=6)
par(mfrow=c(2,1))
par(mar=c(3,3,1,1), mgp=c(1.5,.5,0), tck=-.01)
plot(range(year), range(yreal), type="n", xlab="Year",
      ylab="Day of year (empirical data)", bty="l", main="")
for (j in 1:Nspp){
  lines(year[species==j], yreal[species==j], col="pink3")
 }
plot(range(year), range(y), type="n", xlab="Year", ylab="Day of year (simulated from posterior means)",
     bty="l", main="")
for (j in 1:Nspp)
   lines(year[species==j], y[species==j], col="plum4")
dev.off()

# Create the data using new a and b for each of the species, simshere times
simshere <- 1000
y.sd100.wppint <- matrix(0, ncol=simshere, nrow=Nspp)
for (i in 1:simshere){
    for (n in 1:N){
        s <- species[n]
        ypred[n] <- a[s] + b[s]*year[n] 
    }
  y <- rnorm(N, ypred, sigma_y)
  y.df <- as.data.frame(cbind(y, species))
  y.sd <- aggregate(y.df["y"], y.df["species"], FUN=sd)
  y.sd100.wppint[,i] <- y.sd[,2] 
}

# ... and here's the real data, includes studyid -- which we discussed adding to model
# real.sd <- aggregate(rawlong.nodups["phenovalue"], rawlong.nodups[c("studyid", "spp")], FUN=sd)
real.sd <- aggregate(rawlong.tot2["phenovalue"], rawlong.tot2[c("studyid", "species")],
    FUN=sd)

###############
## Feedbacks ##
###############

if(runmodels){
# See the stan code on this model for notes on what it does
syncmodelhs <- stan("stan/twolevelhierslope.stan", data=c("N","Nspp","y","species","year"),
                   iter=4000, warmup=3000, chains=4, cores=4)
save(syncmodelhs, file="output/syncmodelhs.Rdata")
}

if(!runmodels){
load("output/syncmodelhs.Rdata")
}


# Random slopes only model:
syncmodelhspost <- extract(syncmodelhs)
sigma_bhsmodel <- mean(syncmodelhspost$sigma_b) 
mu_bhsmodel <- mean(syncmodelhspost$mu_b) 

ahsmodel <- colMeans(syncmodelhspost$a) 
bhsmodel <- rnorm(Nspp, mean=mu_bhsmodel, sd=sigma_bhsmodel)


## Redo the prior check above ... 

# extract means for now (other ways to extract the mean)
sigma_y_hsmodel <- mean(syncmodelhspost$sigma_y) 

# Create the data using new a and b for each of the species, simshere times
simshere <- 1000
y.sd100 <- matrix(0, ncol=simshere, nrow=Nspp)
for (i in 1:simshere){
    for (n in 1:N){
        s <- species[n]
        ypred[n] <- ahsmodel[s] + bhsmodel[s]*year[n] 
    }
  y <- rnorm(N, ypred, sigma_y_hsmodel)
  y.df <- as.data.frame(cbind(y, species))
  y.sd <- aggregate(y.df["y"], y.df["species"], FUN=sd)
  y.sd100[,i] <- y.sd[,2] 
}

pdf("graphs/fourpanelforpaper.pdf", height=5.5, width=8)
par(mfrow=c(2,2))
par(mar=c(4,3.5,1,2), mgp=c(1.5,.5,0), tck=-.01)
plot(range(year), range(yreal), type="n", xlab="Year",
      ylab="Day of year \n (empirical data)", bty="l", main="", ylim=c(-10,320))
text(-20, 310, "Empirical data", cex=0.75, col="pink3")
for (j in 1:Nspp){
  lines(year[species==j], yreal[species==j], col="pink3")
 }
plot(range(year), range(y), type="n", xlab="Year", ylab="Day of year \n (simulated from posterior means)",
     bty="l", main="", ylim=c(-10,320))
text(-20, 310, "Simulated data (1X)", cex=0.75, col="plum4")
for (j in 1:Nspp)
   lines(year[species==j], y[species==j], col="plum4")
hist(colMeans(y.sd100.wppint), col="plum3", breaks=20, xlim=c(10,14), 
    main="",
    xlab="Mean SD of response (first model)") # from 1000 sim. datasets (purple) \n versus empirical data (pink)
abline(v = mean(real.sd$phenovalue), col = "pink3", lwd = 2)
text(10.4, 170, "Empircal \n data", cex=0.75, col="pink3")
text(12.15, 160, "Simulated data \n (1000X)", cex=0.75, col="plum4")
hist(colMeans(y.sd100), col="plum3", breaks=20, xlim=c(10,14), 
    main="",
    xlab="Mean SD of response (updated model)")
abline(v = mean(real.sd$phenovalue), col = "pink3", lwd = 2)
text(10.4, 170, "Empircal \n data", cex=0.75, col="pink3")
text(11.8, 160, "Simulated data \n (1000X)", cex=0.75, col="plum4")
dev.off()
