## Started 5 November 2023 ##
## By Lizzie so far ##
## Working on code to source in the supp of the Bayesian workflow paper ##

## housekeeping
rm(list=ls()) 
options(stringsAsFactors = FALSE)

## set working directory if you need to
setwd("/Users/lizzie/Documents/git/projects/misc/miscmisc/bayesianflows/examples/synchrony")

## libraries
library(lme4)

d <- read.csv("output/rawlong.tot2.csv")

modelwanted <- lmer(phenovalue~(year|species), data=d) # singular
modelconverged <- lmer(phenovalue~year+(1|species), data=d)
# summary(modelwanted)