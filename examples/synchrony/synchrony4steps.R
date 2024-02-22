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
library(truncnorm)
library(rstan)
options(mc.cores = parallel::detectCores())

## Step 1: Develop your model (done before we get here)

#################################################
## Simulate data to test code (part of Step 2) ##
#################################################

# Create the species-level parameters
Nspp <- 100 # 88 in the data
mu_doy <- 125
sigma_doy <- 20
mu_shift <- 0.5
sigma_shift <- 1
species_doy <- rnorm(Nspp, mu_doy, sigma_doy)
species_trend <- rnorm(Nspp, mu_shift, sigma_shift)

# Create the overall `error'
sigma_y <- 5

# Keep the parameters together to compare to model output
paramsgiven <- c(mu_doy, mu_shift, sigma_shift, sigma_doy, sigma_y)


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
par(mar=c(3,3,1,1), mgp=c(1.5,.5,0), tck=-.01)
plot(range(year), range(y), type="n", xlab="Year", ylab="Day of year",
     bty="l", main="Test data")
for (sp in 1:Nspp)
  lines(year[species==sp], y[species==sp])


#################################################
## Run the code on test data (part of Step 2) ##
#################################################

fit <- stan("stan/twolevelhierslopeint.stan", data=c("N","y","Nspp","species","year"), iter=1000, chains=4)

# grep stan output
sumer <- summary(fit)$summary
muparams <- sumer[grep("mu", rownames(sumer)), c("25%", "mean", "75%")]
sigmaparams <- sumer[grep("sigma", rownames(sumer)), c("25%", "mean", "75%")]

# compare given versus modeled
paramsgiven
muparams
sigmaparams

spslopes <- sumer[grep("b\\[", rownames(sumer)), "mean"]
plot(spslopes~species_trend, xlab="Given species-level slopes", ylab="Modeled species-level slopes")
abline(0,1)

####################################
## Prior checks (part of Step 2) ##
####################################

# Let's check what the predicted slopes look like
# Iterating over mu and sigma for intercepts and slopes
reps <- 30
mu_doy <- rnorm(reps, 100,30)
sigma_doy <- rtruncnorm(a=0, b=Inf, reps, 0, 20)
mu_shift <- rnorm(reps, 0,5)
sigma_shift <- rtruncnorm(a=0, b=Inf, reps, 0,15)

par(mfrow=c(5,6))
for(i in 1:reps){
    plot(range(year), range(y), xlab="Year", ylab="Day of year",
        xlim=c(-50,40),ylim=c(-50,400), type="n")
    species_doy <- rnorm(Nspp, mu_doy[i], sigma_doy[i])
    species_trend <- rnorm(Nspp, mu_shift[i], sigma_shift[i])
    for(sp in 1:Nspp){
        abline(species_doy[sp], species_trend[sp], col="lightgray")
    }
    abline(mu_doy[i], mu_shift[i], col="black")
}

#######################
## Step 3: real data ##
#######################

# get the data
d <- read.csv("output/rawlong.tot2.csv")

# Formatting for R stan (several ways to do this, this is one)
N <- nrow(d)
y <- d$phenovalue
Nspp <- length(unique(d$species)) #newid is character !
species <- as.numeric(as.factor(d$species))
year <- d$yr1981

if(runmodels){
# See the stan code on this model for notes on what it does
syncmodelhis <- stan("stan/twolevelhierslopeint.stan", data=c("N","Nspp","y","species","year"),
                   iter=4000, warmup=3000, chains=4, cores=4)
save(syncmodelhis, file="output/syncmodelhis.Rdata")
}

if(!runmodels){
load("output/syncmodelhis.Rdata")
}

sumerreal  <- summary(syncmodelhis)$summary
sumerreal[grep("mu", rownames(sumerreal)), c("mean", "2.5%", "25%", "50%", "75%", "97.5%")]
sumerreal[grep("sigma", rownames(sumerreal)), c("mean", "2.5%", "25%", "50%", "75%", "97.5%")]


##########################################
## Posterior predictive checks (Step 4) ##
##########################################

Nreal <- nrow(d)
yreal <- d$phenovalue

# First, plot the real data used in the model
par(mfrow=c(1,2))
par(mar=c(3,3,1,1), mgp=c(1.5,.5,0), tck=-.01)
plot(range(year), range(yreal), type="n", xlab="Year",
     ylab="Day of year", bty="l", main="Raw real data")
for (j in 1:Nspp){
  lines(year[species==j], yreal[species==j])
}
hist(yreal, xlab="Day of year", main="Real data")

# What does a similar plot look like using the model output?
syncmodelhispost <- extract(syncmodelhis) 
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

par(mar=c(3,3,1,1), mgp=c(1.5,.5,0), tck=-.01)
plot(range(year), range(y), type="n", xlab="Year", ylab="Day of year",
    bty="l", main="Data from posterior means")
for (j in 1:Nspp)
  lines(year[species==j], y[species==j])
hist(y, xlab="Day of year", main="Data from posterior means")


# Okay, but that's just one new draw ... PPCs should be done with many draws...
# But then you need to decide on what summary statistics matter because you cannot just look at each plot
# Below I do: SD of y (using the means, should also consider using other draws of the posterior)

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
real.sd <- aggregate(d["phenovalue"], d[c("studyid", "species")],
    FUN=sd)

par(mfrow=c(1,1))
hist(colMeans(y.sd100), col="lightblue", breaks=20, xlim=c(10,14), 
    main="",
    xlab="Mean SD of response from 1000 sim. datasets (light blue) \n versus empirical data (dark blue line)")
abline(v = mean(real.sd$phenovalue), col = "darkblue", lwd = 2)
