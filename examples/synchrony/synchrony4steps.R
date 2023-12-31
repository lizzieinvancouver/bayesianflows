## Started 3 November 2023 ##
## By Lizzie so far ##
## Working on code to source in the supp of the Bayesian workflow paper ##

## Examples in supp of paper
# Figure S2. Estimating phenological change (days/year) using a hinge model and linear model for two interactions. (a) Daphnia spp. (in red) and Perca fluviatillis (in green) over the years 1969 to 2008 and (b) Operophtera brumata (in red) and Parus major (in green) over the years 1961 to 2008. Solid lines represent hinge model and dashed lines represent linear model. The dotted vertical line represents the inflection point of 1981 at year 0.

## housekeeping
rm(list=ls()) 
options(stringsAsFactors = FALSE)

## set working directory if you need to
setwd("/Users/lizzie/Documents/git/projects/misc/miscmisc/bayesianflows/examples/synchrony")

## flags
# You must set to true and RUN the models once for setting this to FALSE to work
runmodels <- FALSE

## libraries
library(rstan)

###########################
## Step 1 simulate data ##
###########################

# Create the species-level parameters
Nspp <- 50
mu_doy <- 125
sigma_doy <- 20
mu_shift <- 0.5
sigma_shift <- 5
species_doy <- rnorm(Nspp, mu_doy, sigma_doy)
species_trend <- rnorm(Nspp, mu_shift, sigma_shift)

# Create the overall `error'
sigma_y <- 5

# Create the data
year_0 <- 1980
n_data_per_species <- round(runif(Nspp, 5, 40))
species <- rep(1:Nspp, n_data_per_species)
N <- length(species)
year <- rep(NA, N)

for (sp in 1:Nspp){
  year[species==sp] <- rev(2009 - 1:(n_data_per_species[sp])) - year_0
}

ypred <- length(N)

for (n in 1:N){
  s <- species[n]
  ypred[n] <- species_doy[s] + species_trend[s]*year[n]
}

y <- rnorm(N, ypred, sigma_y)

# Plot the data

# pdf("testdata.pdf", height=4, width=6)
par(mar=c(3,3,1,1), mgp=c(1.5,.5,0), tck=-.01)
plot(range(year), range(y), type="n", xlab="Year", ylab="Day of year",
     bty="l", main="Test data")
for (sp in 1:Nspp)
  lines(year[species==sp], y[species==sp])
# dev.off()


##############################
## Model fitting for Step 1 ##
##############################

library("rstan")
options(mc.cores = 4)

fit <- stan("stan/twolevelrandomslopeint.stan", data=c("N","y","Nspp","species","year"), iter=1000, chains=4)

# grep stan output
sumer <- summary(fit)$summary
muparams <- sumer[grep("mu", rownames(sumer)), "mean"]
spslopes <- sumer[grep("b\\[", rownames(sumer)), "mean"]

plot(spslopes~species_trend)
abline(0,1)

###########################
## Step 2: Prior checks ##
###########################

# To do ... 

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
syncmodel <- stan("stan/twolevelrandomslope.stan", data=c("N","Nspp","y","species","year"),
                   iter=4000, warmup=3000, chains=4, cores=4)
save(syncmodel, file="output/syncmodel.Rdata")

# See the stan code on this model for notes on what it does
syncmodel2 <- stan("stan/twolevelrandomslopeint.stan", data=c("N","Nspp","y","species","year"),
                   iter=4000, warmup=3000, chains=4, cores=4)
save(syncmodel2, file="output/syncmodel2.Rdata")
}

if(!runmodels){
load("output/syncmodel.Rdata")
load("output/syncmodel2.Rdata")
}

################################
## Posterior predictive checks ##
################################

Nreal <- nrow(rawlong.tot2)
yreal <- rawlong.tot2$phenovalue


# First, plot the real data used in the model
pdf("graphs/realdata_formodel.pdf", height=4, width=6)
par(mar=c(3,3,1,1), mgp=c(1.5,.5,0), tck=-.01)
plot(range(year), range(yreal), type="n", xlab="Year",
     ylab="Day of year", bty="l", main="Raw real data")
for (j in 1:Nspp){
  lines(year[species==j], yreal[species==j])
}
hist(yreal, xlab="Day of year", main="Real data")
dev.off()

# What does a similar plot look like using the model output?
syncmodel2post <- extract(syncmodel2) 
# hist(syncmodel2post$mu_b, xlab="Change per year")

# extract means for now (other ways to extract the mean)
sigma_y <- mean(syncmodel2post$sigma_y) 
sigma_a <- mean(syncmodel2post$sigma_a) 
sigma_b <- mean(syncmodel2post$sigma_b) 
mu_b <- mean(syncmodel2post$mu_b) 
mu_a <- mean(syncmodel2post$mu_a) 

a <- rnorm(Nspp, mean=mu_a, sd=sigma_a)
b <- rnorm(Nspp, mean=mu_b, sd=sigma_b)

N <- Nreal

ypred <- length(N) 
for (n in 1:N){
    s <- species[n]
    ypred[n] <- a[s] + b[s]*year[n]
}
y <- rnorm(N, ypred, sigma_y)

pdf("graphs/onepredictivecheck.pdf", height=4, width=6)
par(mar=c(3,3,1,1), mgp=c(1.5,.5,0), tck=-.01)
plot(range(year), range(y), type="n", xlab="Year", ylab="Day of year",
    bty="l", main="Data from posterior means")
for (j in 1:Nspp)
  lines(year[species==j], y[species==j])
hist(y, xlab="Day of year", main="Data from posterior means")
dev.off()


pdf("graphs/rawvsonepredictivecheck.pdf", height=8, width=6)
par(mfrow=c(2,1))
par(mar=c(3,3,1,1), mgp=c(1.5,.5,0), tck=-.01)
plot(range(year), range(yreal), type="n", xlab="Year",
      ylab="Day of year (empirical data)", bty="l", main="")
for (j in 1:Nspp){
  lines(year[species==j], yreal[species==j])
 }
plot(range(year), range(y), type="n", xlab="Year", ylab="Day of year (simulated from posterior means)",
     bty="l", main="")
for (j in 1:Nspp)
   lines(year[species==j], y[species==j])
dev.off()

# Okay, but that's just one new draw ... PPCs should be done with many draws...
# But then you need to decide on what summary statistics matter because you cannot just look at each plot
# Below I do: SD of y (using the means, I should also consider using other draws of the posterior)

# Create the data using new a and b for each of the species, simshere times
simshere <- 1000
y.sd100 <- matrix(0, ncol=simshere, nrow=Nspp)
for (i in 1:simshere){
    for (n in 1:N){
        s <- species[n]
        ypred[n] <- a[s] + b[s]*year[n] 
    }
  y <- rnorm(N, ypred, sigma_y)
  y.df <- as.data.frame(cbind(y, species))
  y.sd <- aggregate(y.df["y"], y.df["species"], FUN=sd)
  y.sd100[,i] <- y.sd[,2] 
}

# ... and here's the real data, includes studyid -- which we discussed adding to model
# real.sd <- aggregate(rawlong.nodups["phenovalue"], rawlong.nodups[c("studyid", "spp")], FUN=sd)
real.sd <- aggregate(rawlong.tot2["phenovalue"], rawlong.tot2[c("studyid", "species")],
    FUN=sd)

pdf("graphs/retroSDsync.pdf", height=7, width=6)
hist(colMeans(y.sd100), col=rgb(1,0,0,0.5), xlim=c(10,14), 
    main="",
    xlab="Mean(SD of y) from 1000 sim. datasets (pink) versus empirical data (blue)")
abline(v = mean(real.sd$phenovalue), col = "blue", lwd = 2)
dev.off()

# Okay, let's look at other aspects of the model
comppool <- lm(phenovalue~yr1981, data=rawlong.tot2)

# no pooling
uniquespp <- unique(rawlong.tot2$species)
slopefits <- rep(NA< length(uniquespp))
varfits <- rep(NA< length(uniquespp))
intfits <- rep(NA< length(uniquespp))
for(eachsp in 1:length(uniquespp)){
	lmhere <- lm(phenovalue~yr1981, data=subset(rawlong.tot2, species==uniquespp[eachsp]))
	slopefits[eachsp] <- coef(lmhere)[2]
	varfits[eachsp] <- (summary(lmhere)$sigma)**2
	intfits[eachsp] <- coef(lmhere)[1]
}
dr <- data.frame(cbind(uniquespp, slopefits, varfits, intfits))
dr$slopefits <- as.numeric(dr$slopefits)
dr$intfits <- as.numeric(dr$intfits)
dr$varfits <- as.numeric(dr$varfits)

# get 'new' species intercepts and slopes
# this is one way to create fake data from the Stan output to use in the PP check
a <- rnorm(Nspp, mean=mu_a, sd=sigma_a)
b <- rnorm(Nspp, mean=mu_b, sd=sigma_b) 

# compare a few things on this single new dataset
par(mfrow=c(1,2))
hist(b, main="slopes (b) from the stan model with mean from the raw data in blue")
abline(v = mean(dr$slopefits), col = "blue", lwd = 2) # less negative, slopes are generally pooled towards center which makes sense
hist(dr$varfits, main="sigma y (b) from the raw data with sigma_y from the model in blue")
abline(v=sigma_y, col = "blue", lwd = 2) 

par(mfrow=c(1,2))

hist(dr$intfits, breaks=20, main="No pool intercepts", xlab="intercept")
hist(a, breaks=20, main="Partial pool intercepts")


# Random slopes only model:
syncmodelpost <- extract(syncmodel)
sigma_brsmodel <- mean(syncmodelpost$sigma_b) 
mu_brsmodel <- mean(syncmodelpost$mu_b) 

arsmodel <- colMeans(syncmodelpost$a) 
brsmodel <- rnorm(Nspp, mean=mu_brsmodel, sd=sigma_brsmodel)

par(mfrow=c(2,3))
hist(dr$intfits, breaks=20, main="No pool intercepts", xlab="intercept")
hist(a, breaks=20, main="Partial pool intercepts", xlab="intercept")
hist(arsmodel, breaks=20, main="No pool intercepts (with PP slopes)", xlab="intercept")

hist(dr$slopefits, breaks=20, main="No pool slopes", xlab="change over time")
hist(b, breaks=20, main="Partial pool slopes (and intercepts)", xlab="change over time")
hist(brsmodel, breaks=20, main="PP slopes (w/ no pool intercepts)", xlab="change over time")