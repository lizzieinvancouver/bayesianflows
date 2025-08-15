## Started 15 August 2025 ##
## By Lizzie, climate data provided by Victor Van der Meersch ##

## For R Markdown: combine simplest form of extract_and_process.R and lassophenology.R ##

# housekeeping
rm(list=ls()) 
options(stringsAsFactors = FALSE)

wd <- '/Users/lizzie/Documents/git/projects/misc/miscmisc/bayesianflows/examples/lasso'

# libraries
library(geosphere) # get daylength
library(glmnet) # for lasso

# get the climate data 
climdat <- read.csv(file.path(wd, "output/climdatwide.csv"))
tmeandaily <- read.csv(file.path(wd, "output/climdattmean.csv"))

# add daylength
tmeandaily$daylength <- NA
for(i in 1:nrow(tmeandaily)){
  tmeandaily$daylength[i] <- daylength(50, tmeandaily$doy[i])
}

# simulate leafout after 150 
fstar <- 150
lodf <- data.frame(year=unique(tmeandaily$year), 
  loday=rep(NA, length(unique(tmeandaily$year))),
  gdd=rep(NA, length(unique(tmeandaily$year))), 
  daylength=rep(NA, length(unique(tmeandaily$year))))

for(yearhere in unique(tmeandaily$year)) {
  thisyear <- tmeandaily[which(tmeandaily$year==yearhere),]
  leafoutdate <- min(which(cumsum(thisyear[["gddtemp"]]) > fstar))
  lodf$loday[which(lodf$year==yearhere)] <- leafoutdate
  lodf$gdd[which(lodf$year==yearhere)] <- cumsum(thisyear[["gddtemp"]])[leafoutdate]
  lodf$daylength[which(lodf$year==yearhere)] <- thisyear$daylength[leafoutdate]
}

# merge summaries and simulated leafout and daylength
simdat <- merge(climdat, lodf, by="year", all.x=TRUE)


## Now fit lasso regression with all potential variables 
# Create matrix of predictors (X) and response (y)
X <- as.matrix(simdat[, c("tminwinter",  
	"gddspring",  
	"tmeanspring",
	"precspring", 
	"totalprec", 
	"chillwinter",
	"daylength")])
y <- simdat$loday

# Run cross-validated to get lambda and plot results
cv_lasso <- cv.glmnet(X, y, alpha = 1, standardize = TRUE, nfolds = 10)
plot(cv_lasso)

# Find the lambda with smallest cross-validated error and check out the coefficients
best_lambda <- cv_lasso$lambda.min
best_lambda
coef(cv_lasso, s = "lambda.min")

##
## Try again without daylength
X <- as.matrix(simdat[, c("tminwinter",  
	"gddspring",  
	"tmeanspring",
	"precspring", 
	"totalprec", 
	"chillwinter")])
y <- simdat$loday

# Run cross-validated to get lambda and plot results
cv_lasso <- cv.glmnet(X, y, alpha = 1, standardize = TRUE, nfolds = 10)
plot(cv_lasso)

# Find the lambda with smallest cross-validated error and check out the coefficients
best_lambda <- cv_lasso$lambda.min
best_lambda
coef(cv_lasso, s = "lambda.min")

par(mfrow=c(1,2), mgp=c(2, 0.5, 0), tck=-0.01)
plot(loday~tmeanspring, simdat, bty="l", ylab="Leafout day of year", xlab="Mean spring temperature")
plot(loday~daylength, simdat, bty="l", ylab="Leafout day of year", xlab="Daylength")

