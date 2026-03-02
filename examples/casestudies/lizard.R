### Started 15 August 2025 ### 
## By Lizzie ##

## Do we need a super simple example? ##

# housekeeping
rm(list=ls()) 
options(stringsAsFactors = FALSE)

set.seed(7799)

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

# Make overlay histograms of 5 vs 10 reps

reps <- 10 
x2 <- c(rep(0, reps), rep(1, reps))
y2 <- basetoesize + stormeffect*x2 + rnorm(length(x2), 0, sigmay)

hist(y[1:5], breaks=3, xlim=c(2,8), col=rgb(1,0,0,0.5), xlab="Toe size", 
     ylab="Frequency", main="Lizards on islands" , ylim=c(0,6))
hist(y[6:10], breaks=3, xlim=c(0,10), col=rgb(0,0,1,0.5), add=TRUE)

hist(y2[1:10], breaks=3, xlim=c(2,8), col=rgb(1,0,0,0.5), xlab="Toe size", 
     ylab="Frequency", main="Lizards on islands" , ylim=c(0,6))
hist(y2[11:20], breaks=3, xlim=c(0,10), col=rgb(0,0,1,0.5), add=TRUE)


pdf("/Users/lizzie/Documents/git/projects/misc/miscmisc/bayesianflows/examples/lizardexample/lizardgraph.pdf", height=5, width=4)
par(mfrow=c(1,1), mgp=c(2, 0.5, 0), tck=-0.01)
hist(y[1:5], breaks=3, xlim=c(3,9), col=rgb(1,0,0,0.5), xlab="Toe size", 
     ylab="Frequency", main="Lizards on 10 islands" , ylim=c(0,4.5))
hist(y[6:10], breaks=3, xlim=c(0,10), col=rgb(0,0,1,0.5), add=TRUE)
dev.off()