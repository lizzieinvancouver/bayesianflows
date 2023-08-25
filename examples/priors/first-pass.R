source("src/headers.R")

inv.logit <- function(x) exp(x)/(exp(x)+1)
sim.data <- function(n){
  raw.data <- data.frame(x=rnorm(n))
  raw.data$y <- 6 - raw.data$x + rnorm(n)
  return(raw.data)
}
sim.data.broad <- function(n){
  raw.data <- data.frame(x=rnorm(n))
  raw.data$y <- raw.data$x + rnorm(n)
  raw.data$y <- rbinom(nrow(raw.data), 1, inv.logit(raw.data$x))
  return(raw.data)
}
est.coefs <- function(slp.prior, int.prior, data){
  model <- stan_glm(y ~ x, data=data,
                    prior = normal(location=slp.prior, scale=1, autoscale=FALSE),
                    prior_intercept=normal(location=int.prior, scale=1, autoscale=FALSE)
                    )
  return(coef(model))
}
est.coefs.broad <- function(slp.prior, int.prior, data){
  model <- stan_glm(y ~ x, data=data, family=binomial(link="logit"),
                    prior = normal(location=0, scale=slp.prior, autoscale=FALSE),
                    prior_intercept=normal(location=int.prior, scale=1, autoscale=FALSE)
                    )
  return(coef(model))
}


raw.data <- lapply(rep(1000,30), sim.data)
data <- expand.grid(n=c(30,50,100,200,1000), slp.prior=c(1,10,20,40,80), est.slp=NA, est.int=NA, rep=1:30)
for(i in seq_len(nrow(data)))
  data[i,c("est.int","est.slp")] <- est.coefs(data$slp.prior[i], 0, raw.data[[data$rep[i]]][seq_len(data$n[i]),])

raw.broad.data <- lapply(rep(1000,30), sim.data.broad)
broad.data <- expand.grid(n=c(30,50,100,200,1000), slp.prior=c(sqrt(2),3,10,100), est.slp=NA, est.int=NA, rep=1:30)
for(i in seq_len(nrow(broad.data)))
  broad.data[i,c("est.int","est.slp")] <- est.coefs.broad(broad.data$slp.prior[i], 0, raw.broad.data[[broad.data$rep[i]]][seq_len(broad.data$n[i]),])


with(data, plot(abs(slp.prior) ~ n, log="xy", pch=20, col=ifelse(abs(est.slp) >= 2, "red", "black")))


data$thumb <- data$n/data$slp.prior
with(data, plot(est.slp ~ thumb, pch=20, log="x", axes=FALSE, xlab=expression(n / prior), ylab="Estimated slope", cex.lab=1.5, type="n"))
axis(1, at=c(.5,1,10,100,1000), labels=c(.5,1,10,100,1000))
axis(2)
abline(h=-1, col="red", lwd=3)
with(data, points(est.slp ~ thumb, pch=20))


simple <- data.frame(est=as.numeric(with(data, tapply(est.slp, list(n, slp.prior), median))))
simple$prior <- rep(c(1,10,20,40,80), each=5)
simple$n <- rep(c(30,50,100,200,1000), 5)

with(simple, plot())

with(data, table(thumb>10, est.slp<.9))

save.image("wip.RData")








lik.diff <- function(x, model.data){
  data.loglik <- sum(dnorm(model.data, 0, 1, log=TRUE))
  prior.loglik <- dnorm(0, x, 1, log=TRUE)
  return(abs(prior.loglik - data.loglik))
}
optim(0, lik.diff, model.data=data$x, method="Brent", lower=-10000, upper=0)$par



sum(dnorm(data$x, 0, 1, log=TRUE))

dnorm(0, -10, 1, log=TRUE)
