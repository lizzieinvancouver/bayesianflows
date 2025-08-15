### Started 15 August 2025 ### 
## By Lizzie ##

## Do we need a super simple example? ##

# housekeeping
rm(list=ls()) 
options(stringsAsFactors = FALSE)

set.seed(7777)

basetoesize <- 5
stormeffect <- -0.5
sigmay <- 1
reps <- 5
x <- c(rep(0, reps), rep(1, reps))

y <- basetoesize + stormeffect*x + rnorm(length(x), 0, sigmay)

summary(lm(y~x))

df <- data.frame(intercept=numeric(), slope=numeric())
for(i in c(1:100)){
	y <- basetoesize + stormeffect*x + rnorm(length(x), 0, sigmay)
	df[i,] <- coef(lm(y~x))
}

par(mfrow=c(1,2), mgp=c(2, 0.5, 0), tck=-0.01)
hist(df$intercept, bty="l", ylab="Fequency", xlab="Intercept", 
	main="100 simulated datasets")
abline(v=basetoesize, col="dodgerblue", lwd=2)
hist(df$slope, bty="l", ylab="Fequency", xlab="Slope", main="")
abline(v=stormeffect, col="dodgerblue", lwd=2)

nrow(subset(df, slope>0.5))