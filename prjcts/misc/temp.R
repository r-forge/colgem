require(rcolgem)
source('rcolgem.R')
parms_truth <- list(gamma0 = 1
 , gamma1 = 1/7
 , gamma2 = 1/2
 , mu = 1/30
 , b=.036
 , beta0 = 12./10
 , beta1=3./100
 , beta2=9./100
 , S_0=3000
 , I0_0=1, I1_0=0.01, I2_0=0.01
 , m=3, mm=1) 

INFECTEDNAMES <- c('I0', 'I1', 'I2')

births <- rbind( 
  c('parms$beta0 * S * I0 / (S + I0 + I1 + I2)', '0', '0'), 
  c('parms$beta1 * S * I1 / (S + I0 + I1 + I2)', '0', '0'), 
  c('parms$beta2 * S * I2 / (S + I0 + I1 + I2)', '0', '0')
)
rownames(births)=colnames(births)<- INFECTEDNAMES
migrations <- rbind(
  c('0', 'parms$gamma0 * I0', '0'),
  c('0', '0', 'parms$gamma1 * I1'),
  c('0', '0', '0')
)
rownames(migrations)=colnames(migrations) <- INFECTEDNAMES
deaths <- c(
 'parms$mu*I0'
 , 'parms$mu*I1'
 , 'parms$mu*I2 + parms$gamma2 * I2'
)
names(deaths) <- INFECTEDNAMES
nonDemeDynamics <- paste(sep='',
'-parms$mu*S + parms$mu*(S + I0 + I1 + I2)'
 ,'- S * (parms$beta0*I0+parms$beta1*I1+parms$beta2*I2) / (S + I0 + I1 + I2)'
)
names(nonDemeDynamics) <- 'S'


ox <- solve.model(t0=0, t1=50, x0=c(I0=1, I1=.01, I2=.01, S = parms_truth$S_0), births,  deaths, nonDemeDynamics, parms_truth, migrations=migrations, integrationMethod = 'rk4')
plot(ox)


if (F)
{
n <- 300
p <- ox[nrow(ox),INFECTEDNAMES]; p  <- p / sum(p)
sampleTimes <- seq(40, 50, length.out=n) 
sampleStates <- c( rep(1, round(n*p[1])) , rep(2, round(p[2]*n)), rep(3, round(p[3]*n)))
sampleStates <- sample(sampleStates, n , replace=FALSE) #shuffle elements
names(sampleTimes) <- as.character( -n:-1 )
names(sampleStates) <- names(sampleTimes)
print(system.time( bdt <- simulate.binary.dated.tree(births, deaths, nonDemeDynamics
  ,  t0=0
  , x0=c(I0=1, I1=.01, I2=.01, S = parms_truth$S_0)
  , sampleTimes=sampleTimes
  , sampleStates=sampleStates
  ,  migrations=migrations
  ,  parms=parms_truth
  , fgyResolution = 2000
  , integrationMethod = 'rk4')
))
}else
{
	tree <- read.tree(
	   system.file('extdata/hivSimulation.nwk', package='rcolgem'))
	#~ the sample times are the same, because it is a homochronous sample at 50 years
	sampleTimes <- rep(50, length(tree$tip.label))
	names(sampleTimes) <- tree$tip.label
	# create a tree with dated tips and internal nodes, 
	# will infer the sample states from tip labels
	bdt <- binaryDatedTree(tree
	  , sampleTimes
	  , sampleStatesAnnotations=INFECTEDNAMES)
}


print(system.time( print(
	coalescent.log.likelihood( bdt
	 , births,  deaths, nonDemeDynamics
	 , t0 = 0
	 , x0=c(I0=1, I1=.01, I2=.01, S = parms_truth$S_0)
	 , migrations = migrations
	 , parms=parms_truth
	 , fgyResolution=1000
	 , integrationMethod='euler'
	)
)))


if (F)
{ #stable version
source('../../../rcolgem-rforge2/pkg/R/rcolgem.R')
print(system.time( print(
	coalescent.log.likelihood( bdt
	 , births,  deaths, nonDemeDynamics
	 , t0 = 0
	 , x0=c(I0=1, I1=.01, I2=.01, S = parms_truth$S_0)
	 , migrations = migrations
	 , parms=parms_truth
	 , fgyResolution=1000
	 , integrationMethod='euler'
	)
)))
}



if (T)
{library(bbmle)
	obj_fun <- function(lnbeta0, lnbeta1, lnbeta2, t0)
	{
		parms <- parms_truth
		parms$beta0 <- exp(lnbeta0)
		parms$beta1 <- exp(lnbeta1)
		parms$beta2 <- exp(lnbeta2)
		mll <- -coalescent.log.likelihood( bdt
			 , births, deaths, nonDemeDynamics
			 , t0 = t0
			 , x0=c(I0=1, I1=.01, I2=.01, S = parms$S_0)
			 , migrations = migrations
			 , parms=parms
			 , fgyResolution = 1000
			 , integrationMethod = 'rk4'
		)
		# track progress: 
		print(c(mll, exp(c(lnbeta0, lnbeta1, lnbeta2) ), t0) )
		mll
	}

	fit <- mle2(obj_fun
	  , start=list(lnbeta0=log(1), lnbeta1=log(.1), lnbeta2=log(.2), t0=0)
	  , method='Nelder-Mead', optimizer='optim' 
	  ,  control=list(trace=6, reltol=1e-8))
	
}

if (F)
{
	AIC(fit)
	logLik(fit)
	coef(fit)
	exp(coef(fit))
	
	load( system.file('extdata/hivModel0-fit.RData', package='rcolgem') )
	AIC(fit)
	logLik(fit)
	coef(fit)
	exp(coef(fit))
#~ > exp(coef(fit))
#~    lnbeta0    lnbeta1    lnbeta2         t0 
#~ 1.08414135 0.04003263 0.17678659 1.64847045 

#~ > exp(coef(fit))
#~    lnbeta0    lnbeta1    lnbeta2         t0 
#~ 1.19396395 0.03383565 0.05161617 0.15978107 
}

