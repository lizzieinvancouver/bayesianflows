### Started 8 August 2025 ### 
## By Lizzie ###

# housekeeping
rm(list=ls()) 
options(stringsAsFactors = FALSE)

# setwd 
setwd("/Users/lizzie/Documents/git/projects/misc/miscmisc/bayesianflows/examples/lasso")

# libraries
library(glmnet) # for lasso

simdat <- read.csv("output/simuldateddat.csv")
simdat <- simdat[-1,]

###
### Attention! Lasso code from chatGPT (8 August 2025)
###

# Step 1: Create matrix of predictors (X) and response (y)
# Note: glmnet requires X to be a matrix and y to be a vector
X <- as.matrix(simdat[, c("tminwinter",  
	"gddsummer",  
	"precspring", 
	"totalprec", 
	"chillwinter",
	"gdd", 
	"daylength")])
y <- simdat$loday

# Step 2: Run cross-validated Lasso
cv_lasso <- cv.glmnet(X, y, alpha = 1, standardize = TRUE, nfolds = 10)

# Step 3: Plot cross-validation results
plot(cv_lasso)

# Step 4: Find the lambda with minimum cross-validated error
best_lambda <- cv_lasso$lambda.min
best_lambda

# Step 5: Coefficients at best lambda
coef(cv_lasso, s = "lambda.min")

# Step 6: Make predictions
predictions <- predict(cv_lasso, newx = X, s = "lambda.min")

# Step 7: Evaluate performance
# You can use RMSE as a simple metric
rmse <- sqrt(mean((predictions - y)^2))
rmse



## OLDER code
# Get climate data and subset to one site
clim <- read.csv("..//..//..//..//..//treegarden/decsens/analyses/pep_analyses/output/betpen_dailytempsandlo_1950to2010.csv")
head(clim)
table(clim$lat)
clim1 <- subset(clim, lat==50)

# add daylength
clim1$daylength <- NA
clim1$doy <- as.numeric(format(as.Date(clim1$Date, format="%Y-%m-%d") , "%j"))

for(i in 1:nrow(clim1)){
	clim1$daylength[i] <- daylength(50, clim1$doy[i])
}

clim1$gddtemp <- ifelse(clim1$tmean>0, clim1$tmean, 0)
# simulate leafout 
fstar <- 150
lodf <- data.frame(year=unique(clim1$year), 
	loday=rep(NA, length(unique(clim1$year))),
	gdd=rep(NA, length(unique(clim1$year))), 
	daylength=rep(NA, length(unique(clim1$year))))

for(yearhere in unique(clim1$year)) {
	thisyear <- clim1[which(clim1$year==yearhere),]
	leafoutdate <- min(which(cumsum(thisyear[["gddtemp"]]) > fstar))
	lodf$loday[which(lodf$year==yearhere)] <- leafoutdate
	lodf$gdd[which(lodf$year==yearhere)] <- cumsum(thisyear[["gddtemp"]])[leafoutdate]
	lodf$daylength[which(lodf$year==yearhere)] <- thisyear$daylength[leafoutdate]
}

lengthhere <- nrow(lodf)
lodf$var1 <- rnorm(lengthhere, 9, 10)       
lodf$var2 <- rnorm(lengthhere, 2, 5)       
lodf$var3 <- rnorm(lengthhere, 4, 15)       
lodf$var4 <- rnorm(lengthhere, 9, 7)       
lodf$var5 <- rnorm(lengthhere, 9, 1)       
lodf$var6 <- rnorm(lengthhere, 99, 10)   
lodf$var7 <- rnorm(lengthhere, 76, 10)       
lodf$var8 <- rnorm(lengthhere, 0.1, 0.5)       



###
### Attention! Lasso code from chatGPT (8 August 2025)
###

# Step 1: Create matrix of predictors (X) and response (y)
# Note: glmnet requires X to be a matrix and y to be a vector
X <- as.matrix(lodf[, c("var1", "var2", "var3", "var4", "var5",
                      "var6", "var7", "var8",  
                      "gdd", "daylength")])
y <- lodf$loday

# Step 2: Run cross-validated Lasso
cv_lasso <- cv.glmnet(X, y, alpha = 1, standardize = TRUE, nfolds = 10)

# Step 3: Plot cross-validation results
plot(cv_lasso)

# Step 4: Find the lambda with minimum cross-validated error
best_lambda <- cv_lasso$lambda.min
best_lambda

# Step 5: Coefficients at best lambda
coef(cv_lasso, s = "lambda.min")

# Step 6: Make predictions
predictions <- predict(cv_lasso, newx = X, s = "lambda.min")

# Step 7: Evaluate performance
# You can use RMSE as a simple metric
rmse <- sqrt(mean((predictions - y)^2))
rmse


## Try again without GDD
X <- as.matrix(lodf[, c("var1", "var2", "var3", "var4", "var5",
                      "var6", "var7", "var8",  
                       "daylength")])
# Step 2: Run cross-validated Lasso
cv_lasso <- cv.glmnet(X, y, alpha = 1, standardize = TRUE, nfolds = 10)

# Step 3: Plot cross-validation results
plot(cv_lasso)

# Step 4: Find the lambda with minimum cross-validated error
best_lambda <- cv_lasso$lambda.min
best_lambda

# Step 5: Coefficients at best lambda
coef(cv_lasso, s = "lambda.min")

# Step 6: Make predictions
predictions <- predict(cv_lasso, newx = X, s = "lambda.min")

# Step 7: Evaluate performance
# You can use RMSE as a simple metric
rmse <- sqrt(mean((predictions - y)^2))
rmse

##
## Try again without daylength
X <- as.matrix(lodf[, c("var1", "var2", "var3", "var4", "var5",
                      "var6", "var7", "var8",  
                       "gdd")])
# Step 2: Run cross-validated Lasso
cv_lasso <- cv.glmnet(X, y, alpha = 1, standardize = TRUE, nfolds = 10)

# Step 3: Plot cross-validation results
plot(cv_lasso)

# Step 4: Find the lambda with minimum cross-validated error
best_lambda <- cv_lasso$lambda.min
best_lambda

# Step 5: Coefficients at best lambda
coef(cv_lasso, s = "lambda.min")

# Step 6: Make predictions
predictions <- predict(cv_lasso, newx = X, s = "lambda.min")

# Step 7: Evaluate performance
# You can use RMSE as a simple metric
rmse <- sqrt(mean((predictions - y)^2))
rmse