#~ SDE version of model in clusterSizeVignette.R
setwd('/Users/Oliver/workspace_sandbox/rcolgem')
source('rcolgem.R')

#~ parms_truth <<- list(gamma0 = 1, gamma1 = 1/7, gamma2 = 1/2, mu = 1/30, b=.036, beta0 = 0.775, beta1=0.08, beta2=0.08, S0=2500, alpha = .05) 
parms_truth <<- list(gamma0 = 1, gamma1 = 1/7, gamma2 = 1/2, mu = 1/30, b=.036, beta0 = 1+1/30, beta1=(1+1/30)/10, beta2=(1+1/30)/10, S_0=5000, I0_0=1, I1_0=1, I2_0=1, alpha = 4) 
FGYPARMS <<- parms_truth
INFECTEDNAMES <<- c('I0', 'I1', 'I2')

#~ the 'skeleton' functions are deterministic functions of system state
F.skeleton <<- function(t, x = NA)
{
if (is.na(x[1]) )  x <- x.interp(t) 
N <- sum(x) 
I <- sum(x[INFECTEDNAMES])
lambda <- exp(-FGYPARMS$alpha * I/N) * c( FGYPARMS$beta0*x['I0'], FGYPARMS$beta1*x['I1'] , FGYPARMS$beta2*x['I2'])  * x['S']/N
unname( cbind( lambda, matrix(0, nrow=3, ncol=2) ) )
}


G. <<- function(t, x = NA)
{
if (is.na(x[1]) )  x <- x.interp(t) 
g <- matrix(0, nrow=3, ncol=3)
g[1,2] <- FGYPARMS$gamma0 * x['I0']
g[2,3] <- FGYPARMS$gamma1 * x['I1']
g
}


#~ solve using discrete euler method
timestep = .1
times <- seq(0, 50, by=timestep)

step.x <- function(t, x, parms, ...) 
{
with(parms, {
lambda <-  F.skeleton(t, x = x)
lambdaTimestep <- lambda * timestep
sumlambda <- sum(lambda) 
sumlambdaTimestep <- sum(lambdaTimestep) 
rsumlambda <- abs( sumlambda*timestep + sqrt(sumlambda) * rnorm(1, mean=0, sd=sqrt(timestep)) )
lambda <- lambdaTimestep * rsumlambda/sumlambdaTimestep
N <- sum(x)
list( x = x + c(
		S = -rsumlambda + timestep*b * N - timestep * mu * x['S'], 
		I0 = rsumlambda - timestep * (mu+gamma0)*x['I0'],
		I1 = timestep*gamma0*x['I0'] - timestep*(mu+gamma1)*x['I1'],
		I2 = timestep*gamma1*x['I1'] - timestep*(mu+gamma2)*x['I2']) 
	, lambda = lambda )
})
}


solve.model.set.fgy <- function(parameters, shouldPlot=FALSE){ 
	FGYPARMS <<- parameters
	lambdas <- list()
	x <-  c(S = parameters$S_0, I0=parameters$I0_0, I1=parameters$I1_0, I2=parameters$I2_0)
	o <- c(0, x)
	for (itimes in 2:length(times)){
		t <- times[itimes]
		x_lambda <- step.x( t, x, parameters)
		x <- x_lambda[['x']]
		lambdas[[itimes]] <- x_lambda[['lambda']]
		o<- rbind(o, c(t, x))
	}
	x.interps<-  sapply( names(x), function(n) approxfun(o[,1], o[,n], yleft=o[1,n], yright=o[nrow(o),n]) )
	x.interp <<- function(t) { sapply(x.interps, function(interp) interp(t) ) }  
	Y. <<- function(t) x.interp(t)[INFECTEDNAMES]
	F.interps <-  list()
	for (k in 1:m){				#TODO not defined 
		F.interps[[k]] <- list()
		for (l in 1:m){
			lambdaskl <- sapply(2:length(times), function(itime) lambdas[[itime]][k,l])
			F.interps[[k]][[l]] <- approxfun( times[1:(length(times)-1)], lambdaskl/timestep, method='constant', rule=2)
		}
	}
	F. <<- function(t.) t( sapply(1:m, function(k)  sapply(1:m, function(l) F.interps[[k]][[l]](t.)) ) )
	if (shouldPlot){
		class(o) <- 'deSolve'
		#X11()
		pdf( 'clusterSize_stochasticModel_vignette/epidemic.pdf')
		plot(o)
		dev.off()
	}
}

# simulate coalescent tree with true parameters 
solve.model.set.fgy(parms_truth, shouldPlot=TRUE) 
phi = .50 # sample fraction
sampleTime <- 50
Y.sampleTime <- Y.(sampleTime)
stateIndices <- rep( 1:3, round( phi * Y.sampleTime ) ) # sample each of three states in proportion to size
n <- length(stateIndices)
sampleTimes <- rep(sampleTime, n)
sampleStates <- diag(3)[stateIndices,]
coalescentTree_time <- system.time({
  bdt <- simulatedBinaryDatedTree(sampleTime, sampleStates, discretizeRates=TRUE) 
})

# calculate empirical stats
heights <- seq(0, 50, length.out=50)
empiricalMoments = eM <- calculate.cluster.size.moments.from.tree(bdt, heights)


# comparison of model and empirical trees
comparison.plots <- function(eM, mM)
{
	X11()
	par(mfrow=c(3,3))
	for (i in 1:3){
		for (j in 1:3){
			plot( mM$heights, mM$M[i,j,], type='lines', log='y', main=paste(i,j), 
			  ylim=c(1,
					 max(c(mM$M[i,j,], eM[i,j,]) , na.rm = TRUE) )
			  , col='red')
			lines( heights, eM[i,j,])
		}
	}
	X11() 
	plot( -mM$heights, rowSums( mM$A ), type='l' , col='red', main='Lineages through time')
	ltt.lines(bdt) 
}
comparison.plots2 <- function(eM, nsims,parms=parms_truth)
{
	mM <- list()
	for (isim in 1:nsims){
		mM_time <- system.time({  modelMoments = mM_i <- calculate.cluster.size.moments.from.model(sampleTime, sampleStates , timeResolution = 50, discretizeRates=TRUE, fgyResolution = 100 , integrationMethod = 'rk4')
		})
	}
	mM <- lapply(1:nsims, function(isim) {
	  solve.model.set.fgy(parms, shouldPlot=FALSE) ;
	  calculate.cluster.size.moments.from.model(sampleTime, sampleStates , timeResolution = 50, discretizeRates=TRUE, fgyResolution = 100 , integrationMethod = 'rk4');
	  }
	 )
	X11()
	par(mfrow=c(4,3))
	for (i in 1:3){
		for (j in 1:3){
			plot( heights, eM[i,j,], type='lines', log='y', main=paste(i,j), 
			  ylim=c(1,
					 max(c(mM$M[i,j,], eM[i,j,]) , na.rm = TRUE) )
			  , col='black', 
			  xlab='height')
			for (isim in 1:nsims){
				lines(mM[[isim]]$heights, mM[[isim]]$M[i,j,], col='red')
			}
		}
	}
	
	ltt.plot(bdt) 
	for (isim in 1:nsims){ 
		lines( -mM[[isim]]$heights, rowSums( mM[[isim]]$A ), col='red')
	}
}
comparison.plots3<- function(eM, nsims,parms=parms_truth)
{ # deltas
	mM <- list()
	for (isim in 1:nsims){
		mM_time <- system.time({  modelMoments = mM_i <- calculate.cluster.size.moments.from.model(sampleTime, sampleStates , timeResolution = 50, discretizeRates=TRUE, fgyResolution = 100 , integrationMethod = 'rk4')
		})
	}
	mM <- lapply(1:nsims, function(isim) {
	  solve.model.set.fgy(parms, shouldPlot=FALSE) ;
	  calculate.cluster.size.moments.from.model(sampleTime, sampleStates , timeResolution = 50, discretizeRates=TRUE, fgyResolution = 100 , integrationMethod = 'rk4');
	  }
	 )
	#X11()
	pdf( 'clusterSize_stochasticModel_vignette/delta_Mij.pdf')
	par(mfrow=c(3,3))
	nheights <- length(heights)
	for (i in 1:3){
		for (j in 1:3){
			emdeltas <- eM[i,j,2:nheights]-eM[i,j,1:(nheights-1)]
			plot( heights[2:nheights], emdeltas, type='points', log='y', main=paste(i,j), 
			  ylim=c(.01,
					 max(emdeltas +10 , na.rm = TRUE) )
			  , col='black', 
			  xlab='height')
			for (isim in 1:nsims){
				mmdeltas <- mM[[isim]]$M[i,j,2:nheights] - mM[[isim]]$M[i,j,1:(nheights-1)]
				lines(mM[[isim]]$heights[2:nheights], mmdeltas, col='red')
			}
		}
	}
	dev.off()
	
	#residuals
	#X11()
	pdf( 'clusterSize_stochasticModel_vignette/residuals.pdf')
	par(mfrow=c(3,3))
	print('test')
	for (i in 1:3){
		for (j in 1:3){
			emdeltas <- log( eM[i,j,2:nheights]-eM[i,j,1:(nheights-1)] )
			isim <- 1
			mmdeltas <- log(mM[[isim]]$M[i,j,2:nheights] - mM[[isim]]$M[i,j,1:(nheights-1)])
			plot( heights[2:nheights], mmdeltas-emdeltas, type='points', main=paste(i,j), 
			  col='black', 
			  xlab='height')
			for (isim in 2:nsims){
				mmdeltas <- log( mM[[isim]]$M[i,j,2:nheights] - mM[[isim]]$M[i,j,1:(nheights-1)] )
				points(mM[[isim]]$heights[2:nheights], mmdeltas-emdeltas, col='black')
			}
		}
	}
	dev.off()
	
	#histogram of residuals
	#X11()
	pdf( 'clusterSize_stochasticModel_vignette/residuals_hist.pdf')
	par(mfrow=c(3,3))
	for (i in 1:3){
		for (j in 1:3){
			mresiduals <- c()
			emdeltas <- log( eM[i,j,2:nheights]-eM[i,j,1:(nheights-1)] )
			isim <- 1
			mmdeltas <- log(mM[[isim]]$M[i,j,2:nheights] - mM[[isim]]$M[i,j,1:(nheights-1)])
			mresiduals <- c(mmdeltas - emdeltas, mresiduals)
			for (isim in 2:nsims){
				mmdeltas <- log( mM[[isim]]$M[i,j,2:nheights] - mM[[isim]]$M[i,j,1:(nheights-1)] )
				mresiduals <- c(mmdeltas - emdeltas, mresiduals)
			}
			hist(mresiduals,main=paste(i,j) )
		}
	}
	dev.off()
}

comparison.plots3(eM , 20)



# calculate model stats
#~ mM_time <- system.time({
#~   modelMoments = mM <- calculate.cluster.size.moments.from.model(sampleTime, sampleStates , timeResolution = 50, discretizeRates=TRUE, fgyResolution = 100 , integrationMethod = 'rk4')
#~ })
#~ comparison.plots(eM, mM) 

#~ now do a comparison with model parameters different from true parameters 
#~ parameters <- parms_truth
#~ parameters$beta0 <- parms_truth$beta0/2
#~ parameters$beta1 <- parms_truth$beta1*5
#~ solve.model.set.fgy(parameters)
#~ modelMoments = mM <- calculate.cluster.size.moments.from.model(sampleTime, sampleStates , timeResolution = 50, discretizeRates=TRUE, fgyResolution = 100 , integrationMethod = 'rk4')
#~ comparison.plots(eM, mM) 


#~ compare computation times with likelihood
#~ solve.model.set.fgy(parms_truth)
#~ likelihood_time <- system.time({
#~   print( paste( 'likelihood', coalescent.log.likelihood(bdt, integrationMethod = 'rk4', finiteSizeCorrections=TRUE, maxHeight=0, discretizeRates=TRUE, fgyResolution = 100) ))
#~ })
#~ print('coalescent tree time')
#~ print(coalescentTree_time)
#~ print('running time moments')
#~ print(mM_time)
#~ print('running time likelihood')
#~ print(likelihood_time)
#~ [1] "likelihood -2225.51092458556"
#~ [1] "coalescent tree time"
#~    user  system elapsed 
#~   0.704   0.000   0.701 
#~ [1] "running time moments"
#~    user  system elapsed 
#~   0.216   0.000   0.216 
#~ [1] "running time likelihood"
#~    user  system elapsed 
#~  15.121   0.004  19.464 
