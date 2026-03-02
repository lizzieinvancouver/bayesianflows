### Started 8 August 2025 ### 
## By Lizzie ###

# housekeeping
rm(list=ls()) 
options(stringsAsFactors = FALSE)

# setwd 
setwd("/Users/lizzie/Documents/git/projects/misc/miscmisc/bayesianflows/examples/lasso")

# libraries
library(glmnet) # for lasso

# get the data 
simdat <- read.csv("output/simuldateddat.csv")


###
### Attention! Lasso code from chatGPT (8 August 2025)
###

# Step 1: Create matrix of predictors (X) and response (y)
# Note: glmnet requires X to be a matrix and y to be a vector
X <- as.matrix(simdat[, c("tminwinter",  
	"gddspring",  
	"tmeanspring",
	"precspring", 
	"totalprec", 
	"chillwinter",
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

##
## Try again without GDD
X <- as.matrix(simdat[, c("tminwinter",  
	"tmeanspring",
	"precspring", 
	"totalprec", 
	"chillwinter",
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
X <- as.matrix(simdat[, c("tminwinter",  
	"precspring", 
	"tmeanspring",
	"totalprec", 
	"chillwinter",
	"gddspring")])
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