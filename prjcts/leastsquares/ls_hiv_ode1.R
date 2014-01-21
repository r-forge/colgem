require(dfoptim)
require(ape)
require(deSolve)
require(rcolgem)

print(paste( date(), 'start') )

#~ parms_truth <<- list(gamma0 = 1, gamma1 = 1/7, gamma2 = 1/2, mu = 1/30, b=.036, beta0 = 0.775, beta1=0.08, beta2=0.08, S0=2500, alpha = .05) 
parms_truth <<- list(gamma0 = 1, gamma1 = 1/7, gamma2 = 1/2, mu = 1/30, b=.036, beta = 6/10, w1=1/10, w2=3/10, S_0=5000, I0_0=1, I1_0=0.01, I2_0=0.01, alpha = 5, m=3) 
INFECTEDNAMES <- c('I0', 'I1', 'I2')
ESTNAMES <- c('alpha',  'beta', 'w1', 'w2', 'I0_0')


################################################
#
#	script helper functions
#
################################################

	
	
	################################################
	#
	#	define F G solve.model.set.fgy that are assumed in rcolgem.R
	#
	################################################		
	#~ the 'skeleton' functions are deterministic functions of system state
	F.skeleton <- function(t, parms, x = NA,  ...)
	{
		if (is.na(x[1]) )  x <- x.interp(t) 
		N <- sum(x) 
		I <- sum(x[INFECTEDNAMES])
		inc_weights <- c( 1*x['I0'], parms$w1*x['I1'] , parms$w2*x['I2'])  
		lambda <- exp(-parms$alpha * I/N) * (parms$beta * I * x['S']/N) * inc_weights / sum(inc_weights)
		unname( cbind( lambda, matrix(0, nrow=3, ncol=2) ) )
	}
	G. <- function(t, x = NA, parms = NA, ...)
	{
		if (is.na(x[1]) )  x <- x.interp(t) 
		if (is.na(parms[1]) ) parms <- parms_truth
		g <- matrix(0, nrow=3, ncol=3)
		g[1,2] <- parms$gamma0 * x['I0']
		g[2,3] <- parms$gamma1 * x['I1']
		g
	}	
	#~ solve using discrete euler method
	timestep	<- 1
	times 		<- seq(0, 50, by=timestep)	
	step.x 		<- function(t, x, parms, ...) 
	{ 
		with(parms, {
					lambda <-  F.skeleton(t, parms, x = x)
					lambdaTimestep <- lambda * timestep
					sumlambda <- sum(lambda) 
					sumlambdaTimestep <- sum(lambdaTimestep) 
					rsumlambda <- max( sumlambda*timestep + sqrt(sumlambda) * rnorm(1, mean=0, sd=sqrt(timestep)) , 0)
					lambda <- lambdaTimestep * rsumlambda/sumlambdaTimestep
					N <- sum(x)
					#if (is.na(-rsumlambda + timestep*b * N - timestep * mu * x['S'])) browser()
					list( x = x + c(
									S = -rsumlambda + timestep*b * N - timestep * mu * x['S'], 
									I0 = rsumlambda - timestep * (mu+gamma0)*x['I0'],
									I1 = timestep*gamma0*x['I0'] - timestep*(mu+gamma1)*x['I1'],
									I2 = timestep*gamma1*x['I1'] - timestep*(mu+gamma2)*x['I2']) 
							, lambda = lambda )
				})
	}
	solve.sde.model.set.fgy <- function(parms, shouldPlot=FALSE, dir.name=NA)
	{ 
		lambdas <- list()
		x <-  c(S = parms$S_0, I0=parms$I0_0, I1=parms$I1_0, I2=parms$I2_0)
		names(x) <- c('S', INFECTEDNAMES)
		o <- c(0, x)
		for (itimes in 2:length(times)){
			t <- times[itimes]
			x_lambda <- step.x( t, x, parms)
			x <- x_lambda[['x']]
			lambdas[[itimes]] <- x_lambda[['lambda']]
			o<- rbind(o, c(t, x))
		}
		x.interps<<-  sapply( names(x), function(n) approxfun(o[,1], o[,n], yleft=o[1,n], yright=o[nrow(o),n]) )
		x.interp <<- function(t) { sapply(x.interps, function(interp) interp(t) ) }  
		Y. <- function(t) x.interp(t)[INFECTEDNAMES]
		F.interps <<-  list()
		for (k in 1:parms$m){
			F.interps[[k]] <- list()
			for (l in 1:parms$m){
				lambdaskl <- sapply(2:length(times), function(itime) lambdas[[itime]][k,l])
				
				tryCatch( {F.interps[[k]][[l]] <- approxfun( times[1:(length(times)-1)], lambdaskl/timestep, method='constant', rule=2)},
				  error = function(e) browser() )
			}
		}
		F. <- function(t.) t( sapply(1:parms$m, function(k)  sapply(1:parms$m, function(l) F.interps[[k]][[l]](t.)) ) )
		
		FGY <- list( F. = F., G. = G., Y. = Y. )
		
		if (shouldPlot)
		{
			class(o) <- 'deSolve'
			#X11()
			#pdf( paste(dir.name,'epidemic.pdf', sep='/') )
			plot(o)
			#dev.off()
		}
		FGY
	}
	
	dx <- function(t, x, parms, ...) 
	{
		with(parms, {
					lambda 	<- sum( F.skeleton(t, parms, x = x ))
					N 		<- sum(x)
					list(c(
									S = -lambda + b * N - mu * x['S'], 
									I0 = lambda - (mu+gamma0)*x['I0'],
									I1 = gamma0*x['I0'] - (mu+gamma1)*x['I1'],
									I2 = gamma1*x['I1'] - (mu+gamma2)*x['I2']
							)) 
				})
	}	
	#	solve ODE and return solution in Y.
	solve.ode.model.set.fgy <- function(parameters)
	{ 
		times 		<- seq(from=0, to=50, length.out = 100)
		x0 			<- c(S = parameters$S_0, I0=parameters$I0_0, I1=parameters$I1_0, I2=parameters$I2_0) 
		names(x0) <- c('S', INFECTEDNAMES)
		o 			<- ode( x0, times, func=dx, parms=parameters, method='adams')
		x.interps	<<-  sapply( names(x0), function(n) approxfun(o[,1], o[,n], yleft=o[1,n], yright=o[nrow(o),n]) )
		x.interp 	<<- function(t) { sapply(x.interps, function(interp) interp(t) ) }  
		
		Y. 			<- function(t) x.interp(t)[INFECTEDNAMES]
		F. 			<- function(t) F.skeleton(t, parameters, x = x.interp(t))
		
		list( F. = F., G. = G., Y. = Y. )
	}
	
	


################################################
#
#	start script
#
################################################


# simulate coalescent tree with true parameters 
FGY <- solve.sde.model.set.fgy(parms_truth, shouldPlot=FALSE ) 
times0 		<- seq(from=0, to=50, length.out = 100)
Ys <- t( sapply(times0, FGY$Y.) )
phi				<- .750 # sample fraction
sampleTime 		<- 50
Y.sampleTime 	<- FGY$Y.(sampleTime)
m				<- 3
stateIndices 	<- rep( 1:m, round( phi * Y.sampleTime ) ) # sample each of three states in proportion to size
sampleSize 		<- length(stateIndices)
sampleTimes 	<- rep(sampleTime, sampleSize)
sampleStates 	<- diag(m)[stateIndices,]
treeSimTime <- system.time( {bdt <- simulatedBinaryDatedTree(sampleTime, sampleStates, FGY = FGY, discretizeRates=TRUE) })
print(paste(date(), 'simulated tree'))
print(treeSimTime)


# calculate empirical stats
heights <- c( seq(0,4,by=.5), seq(6, 40, by=2)) #TODO how to pick range & resolution is an open question
nheights <- length(heights)
empClusterCalcTime<- system.time( {.eM <- calculate.cluster.size.moments.from.tree(bdt, heights) })
empiricalMoments= eM <- log( 1+ .eM)
eAi_heights <-  log( 1 + sapply( heights, function(h) length(.extant.at.height(h, bdt)) ) )
print(eAi_heights)
print(paste(date(), 'calculated tree moments'))
print(empClusterCalcTime)

#~ calculate Deltas
eMij_heights <- c()
emdeltas <- c()
eadeltas <- eAi_heights[2:nheights] - eAi_heights[1:(nheights-1)]
for (i in 1:3){
	for (j in i:3){
		emdeltas <- cbind( emdeltas,  eM[i,j,2:nheights]-eM[i,j,1:(nheights-1)] )
		eMij_heights <- cbind( eMij_heights, eM[i,j,] )
	}
}
#~ eMij_heights_max <- sapply(1:ncol(eMij_heights), function(i) max(eMij_heights[,i]) )
#~ these are the data to be fitted:
#~ X <- cbind( emdeltas, eadeltas)
X <- cbind( eMij_heights,  eAi_heights ) 


#~ define moment condition
#~ should estimate alpha, beta0, beta1, beta2, and I0; work with log-transformed parameter values
.gobj <- function(theta, x, timeInfo=FALSE, ALTOUT=FALSE)
{
	parms <- parms_truth
	parms[ESTNAMES] <- exp(theta) 
	FGY <- solve.ode.model.set.fgy(parms)
	if (timeInfo)
	{
		modelMomentsTime <- system.time( {
		  .mM <- calculate.cluster.size.moments.from.model(sampleTime, sampleStates , FGY=FGY,
			timeResolution = 200, discretizeRates=TRUE, fgyResolution = 100, integrationMethod = 'euler')
		  })
		print('model moments time')
		print(modelMomentsTime)
	} else
	{
		.mM <- calculate.cluster.size.moments.from.model(sampleTime, sampleStates , FGY=FGY,
			timeResolution = 200, discretizeRates=TRUE, fgyResolution = 100, integrationMethod = 'euler')
	}
	mM <- log( 1 +  .mM$M )
	mA <-   1 + .mM$A  
	#if (sum( is.nan( mM )) > 0) return(matrix(Inf, nrow=nrow(x), ncol=ncol(x)) ) #browser()
	mmdeltas <- c()
	.mMij_heights <- c() 
	.A_heights <- c()
	k <- 0
	for (i in 1:3){
		A_heights <- tryCatch( { A_heights <- approx( .mM$heights, mA[,i], xout=heights, rule=2)$y}, 
		  error = function(e) {rep(0, length(heights)); })
		.A_heights <- cbind(  .A_heights, A_heights )
		for (j in i:3){
			k <- k+1
			mMij_heights <- tryCatch( {approx( .mM$heights, mM[i,j,], xout=heights, rule=2)$y }, 
			  error = function(e) rep(0, length(heights)) )
			.mMij_heights <- cbind( .mMij_heights, mMij_heights )
			mmdeltas <- cbind( mmdeltas, mMij_heights[2:nheights] - mMij_heights[1:(nheights-1)] )
		}
	}
	rsA <- log( 1 + rowSums(.A_heights) )
	if (ALTOUT) return( cbind(.mMij_heights, rsA))
	cbind( .mMij_heights, rsA) - x
}



theta0 <- log( unlist( parms_truth[ESTNAMES] ))  # also set below
gg <- .gobj(theta0, X, timeInfo = TRUE)
print(gg)
nstats <- length(gg)

#~ calculate likelihood for time comparison
#~ lik_time <- system.time( {ll <- coalescent.log.likelihood(bdt, FGY=FGY, integrationMethod = 'euler', finiteSizeCorrections=TRUE, maxHeight=0, discretizeRates=TRUE, fgyResolution = 100)} )
#~ print(lik_time)
#~ print(ll)

if (FALSE)
{ # compare stoch and det trajectories
	FGY_det <- solve.ode.model.set.fgy(parms_truth)
	ydet <- sapply( times0, FGY_det$Y.) 
	plot(times0, ydet[2,], type='lines', col='black')
	for (i in 1:nrow(ydet)) lines(times0, ydet[i,], col='black')
	FGY_stoch <- solve.sde.model.set.fgy(parms_truth)
	ystoch <- sapply( times0, FGY_stoch$Y.) 
	for (i in 1:nrow(ystoch)) lines(times0, ystoch[i,], col='red')
}


if (FALSE)
{ # compare cluster moments at true parameter values
	X11()
	plot( heights, eMij_heights[,1], type='lines', ylim=c(0,8))
	for (i in 2:ncol(eMij_heights)) lines( heights ,eMij_heights[,i] )
	lines( heights , eAi_heights )
	mM <- .gobj(theta0, X, ALTOUT=TRUE)
	for (i in 1:ncol(mM)) { 
		tryCatch( { lines( heights,  mM[,i], lty='dashed' , col='red') }, 
		  error = function(e) browser() )
	}
	
	
	X11()
	plot( heights, eMij_heights[,1], type='lines', ylim=c(0,8))
	for (i in 2:ncol(eMij_heights)) lines( heights ,eMij_heights[,i] )
	
	FGY <- solve.sde.model.set.fgy(parms_truth)
	bdt2 <- simulatedBinaryDatedTree(sampleTime, sampleStates, FGY = FGY, discretizeRates=TRUE) 
	.eM2 <- calculate.cluster.size.moments.from.tree(bdt2, heights) 
	eM2 <- log( 1+ .eM2)
	eMij_heights2 <- c()
	for (i in 1:3){
		for (j in i:3){
			eMij_heights2 <- cbind( eMij_heights2, eM2[i,j,] )
		}
	}
	for (i in 1:ncol(eMij_heights2)) { lines( heights,  eMij_heights2[,i], lty='dashed' , col='red') }
}


sample.initial.conditions <- function()
{
	theta0 <- c( 
	  runif(1, 1, 7)
	  , runif(1, parms_truth[['beta']]/2, parms_truth[['beta']]*2)
	  , runif(2, .05, 1.5)
	  , runif(1, .5, 2) )
	log(theta0)
}

#~ fit the model
if (TRUE)
{
	.obj1 <- function(thet, x, w,  gf)
	{
		# calculates weighted RSS
		gt <- as.vector(gf(thet, x))
		as.vector(  gt%*% (w %*% gt) )
	}
	print(paste(date(), 'fitting model'))
	
	w=diag(nstats) #identity matrix for OLS
	
	fits <- lapply( 1:10, function(i)
	{
		print(paste(date(), 'start fit', i))
		theta0 <- sample.initial.conditions()
		print(paste('initial conditions:', exp(theta0)))
		optim(theta0, .obj1, x = X, w = w, gf = .gobj,  control=list(trace=5,reltol=1e-04,maxit=100))
		#hjk(theta0, .obj1, x = X, w = w, gf = gobj, nsim=1, control=list(trace=5,tol=1e-02,maxfeval=300))
	})
	print(paste(date(), 'done with fits iter 1'))
	bestfit_i <- which.min( sapply(1:10, function(i) fits[[i]]$value) )
	fit <- fits[[bestfit_i]]
	if (fit$convergence!=0)
	{
		fit <- optim(fit$par, .obj1, x = X, w = w, gf = .gobj,  control=list(trace=5,reltol=1e-04,maxit=100))
 		#fit <- nmk(fit$par, .obj1, x = X, w = w, gf = gobj,  control=list(trace=5,tol=1e-02,maxfeval=300))
	}
	
	print(fit)
	results <- rbind( unlist(parms_truth[ESTNAMES])
	  , exp(fit$par))
	colnames(results) <- ESTNAMES
	print(results)
	
	g_truth <- .gobj( log(unlist(parms_truth[ESTNAMES])), X)
	g_fitted <- .gobj( fit$par, X)
	print(paste('RSS truth', sum(g_truth^2)))
	print(paste('RSS fitted', sum(g_fitted^2)))
	
	
	X11(); plot(times0, Ys[, 'I1'], type='lines')#, ylim = c(0, 800))
	for (n in INFECTEDNAMES){ lines(times0, Ys[,n], lty='solid' ) }
	parameters <- parms_truth
	parameters[ESTNAMES] <- exp(fit$par)
	for (i in 1:10){
		FGY <- solve.ode.model.set.fgy(parameters)
		o <- t(sapply(times0, FGY$Y.))
		colnames(o) <- INFECTEDNAMES
		for (n in INFECTEDNAMES) { lines( times0,  o[,n], lty='dashed' ) }
	}
	
#~ 	X11()
#~ 	hist(g_fitted, main='gfitted')
#~ 	
#~ 	X11()
#~ 	hist(g_truth, main='truth')
	
	
	#~ plot fitted versus actual moment(height)
	X11()
	plot( heights, eMij_heights[,1], type='lines', ylim=c(0,10.1))
	for (i in 2:ncol(eMij_heights)) lines( heights ,eMij_heights[,i] )
	lines( heights ,eAi_heights )
	for (k in 1:5){
		mM <- .gobj(fit$par, X, ALTOUT=TRUE)
		for (i in 1:ncol(mM)) { lines( heights,  mM[,i], lty='dashed' , col='red') }
	}
	
	
	
	paf.through.time <- function(theta)
	{
		parameters <- parms_truth
		parameters[ESTNAMES] <- exp(theta)
		FGY <- solve.ode.model.set.fgy(parameters)
		o <- t( sapply( heights, function(h) FGY$F.(50-h)[,1]/sum(FGY$F.(50-h)[,1])  ) )
		colnames(o) <- INFECTEDNAMES
		o
	}
	
	X11()
	pcol=c(I0='red', I1='blue', I2='green')
	true_paf <-  paf.through.time( theta0 )
	plot(heights, true_paf[, 'I0'], type='lines', ylim = c(0, 1))
	for (n in INFECTEDNAMES){ lines(heights, true_paf[,n], lty='solid', col=pcol[n] ) }
	fitted_paf <- paf.through.time( fit$par)
	for (n in INFECTEDNAMES) { lines( heights,  fitted_paf[,n], lty='dashed' , col=pcol[n]) }
	print(true_paf[1,])
	print(fitted_paf[1,])
}


#~ 1.7 sec to calc stats from tree
#~ .5 sec to simulate tree
#~ 2.2 total
#~ .5 sec to analytically solve moments for model
