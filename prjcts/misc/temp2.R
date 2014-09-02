require(rcolgem)
source('rcolgem.R')

library(rcolgem)
tree <- read.tree(system.file( 'extdata/sirModel0.nwk', package='rcolgem'))

parms_truth <- list( beta = .00020002, gamma = 1, S0 = 9999, t0 = 0 )

sampleTimes <- rep(12, 75)
names(sampleTimes) <- tree$tip.label
bdt <- binaryDatedTree( tree, sampleTimes=sampleTimes) 
bdt

INFECTEDNAMES <- c('I0', 'I1', 'I2')

births <- c( I = 'parms$beta * S * I' )
deaths <- c( I = 'parms$gamma * I' )
nonDemeDynamics <- c(S = '-parms$beta * S * I')

x0 <- c(I=1, S= unname(parms_truth$S0) )
t0 <- bdt$maxSampleTime - max(bdt$heights) -1 

print( 
system.time( 
print(
  coalescent.log.likelihood(
    bdt
    , births, deaths, nonDemeDynamics
    , t0, x0
    , parms=parms_truth
    , fgyResolution = 1000
    , integrationMethod = 'rk4')
)))


if (T)
{
library(bbmle)
obj_fun <- function(lnbeta, lnI0)
{
	beta <- exp(lnbeta)
	I0 <- exp(lnI0)
	parms <- parms_truth
	parms$beta <- beta
	x0 <- c(I=unname(I0), S = unname(parms$S0) )
	mll <- -coalescent.log.likelihood(
		bdt
		, births, deaths, nonDemeDynamics
		,  t0, x0
		, parms=parms
		, fgyResolution = 500
		, integrationMethod = 'rk4')
	print(paste(mll, beta, I0))
	mll
}

fit <- mle2(
  obj_fun
  , start=list(lnbeta=log(parms_truth$beta*.75), lnI0=log(1))
  , method='Nelder-Mead'
  , optimizer='optim' 
  , control=list(trace=6, reltol=1e-8)
)


}

if (F)
{
	AIC(fit)
	logLik(fit)
	coef(fit)
	exp(coef(fit))
	# how biased is the estimate? 
	exp(coef(fit)['lnbeta']) - parms_truth$beta

	load( system.file('extdata/sirModel0-fit.RData', package='rcolgem') )
	AIC(fit)
	logLik(fit)
	coef(fit)
	exp(coef(fit))
	# how biased is the estimate? 
	exp(coef(fit)['lnbeta']) - parms_truth$beta

}
