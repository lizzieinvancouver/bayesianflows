## Started Halloween 2025 (Gobabeb) ##
## Working off Will's  misspecified-regression.rmd ##
## ... to make a figure for the paper ##

# Simulate one lineage
apples.x <- rnorm(n=100, mean=-1, sd=1)
apples.y <- rnorm(n=100, mean=-1, sd=1)

# Simulate another lineage
oranges.x <- rnorm(n=100, mean=1, sd=1)
oranges.y <- rnorm(n=100, mean=1, sd=1)

# Merge the data
data <- data.frame(x=c(apples.x,oranges.x), 
  y=c(apples.y,oranges.y), fruit=rep(c("apples","oranges"),each=100))

pdf("/Users/lizzie/Documents/git/projects/misc/miscmisc/bayesianflows/examples/misspecifiednonidentified/nonidentgraph.pdf", 
		height=5, width=4)
par(mgp=c(2, 0.5, 0), tck=-0.01)
with(data, plot(y~x, type="n"))
points(data=subset(data, fruit=="apples"), y~x, col="red")
points(data=subset(data, fruit=="oranges"), y~x, col="orange")
legend("topleft", col=c("red", "orange"), legend=c("clade 1", "clade 2"), pch=1, bty="n")
dev.off()