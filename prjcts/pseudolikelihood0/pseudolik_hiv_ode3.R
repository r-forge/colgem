#~ based on pseudolik...2, this version keeps alpha fixed, uses stats4::mle
require(dfoptim)
require(ape)
require(deSolve)
#~ require(stats4)
require(bbmle)
#~ source('rcolgem.R')
source('rcolgem-devel.R')

print(paste( date(), 'start') )

#~ parms_truth <<- list(gamma0 = 1, gamma1 = 1/7, gamma2 = 1/2, mu = 1/30, b=.036, beta0 = 0.775, beta1=0.08, beta2=0.08, S0=2500, alpha = .05) 
parms_truth <<- list(gamma0 = 1, gamma1 = 1/7, gamma2 = 1/2, mu = 1/30, b=.036, beta0 = 1+1/30, beta1=(1+1/30)/10, beta2=(1+1/30)/10, S_0=5000, I0_0=1, I1_0=0.01, I2_0=0.01, alpha = 4, m=3) 
INFECTEDNAMES <<- c('I0', 'I1', 'I2')


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
	F.skeleton <- function(t,parms , x = NA,  ...)
	{
		if (is.na(x[1]) )  x <- x.interp(t) 
		N <- sum(x) 
		I <- sum(x[INFECTEDNAMES])
		lambda <- exp(-parms$alpha * I/N) * ( x['S']/N) * c( parms$beta0*x['I0'], parms$beta1*x['I1'] , parms$beta2*x['I2'])  
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
					lambda <-  F.skeleton(t, parms,  x = x)
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
					lambda 	<- sum( F.skeleton(t, parms, x = x))
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
phi				<- .50 # sample fraction
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
heights <- seq(0, max(bdt$heights), length.out=21) #TODO how to pick range & resolution is an open question
nheights <- length(heights)

ESTNAMES <- c( 'beta0', 'beta1', 'beta2', 'I0_0')
theta0 <- log( unlist( parms_truth[ESTNAMES] ))  # also set below


#~ calculate likelihood for time comparison
#~ lik_time <- system.time( {ll <- coalescent.log.likelihood(bdt, FGY=FGY, integrationMethod = 'euler', finiteSizeCorrections=TRUE, maxHeight=0, discretizeRates=TRUE, fgyResolution = 100)} )
#~ print(lik_time)
#~ print(ll)


if (FALSE)
{ # compare stoch and det trajectories
	FGY_stoch <- solve.sde.model.set.fgy(parms_truth)
	ystoch <- sapply( times0, FGY_stoch$Y.) 
	plot(times0, colSums(ystoch), type='lines', col='black')
	for (i in 1:nrow(ystoch)) lines(times0, ystoch[i,], col='red')
}


if(FALSE)
{
	etime <- system.time( {heights_entanglements <- calculate.coalescent.entanglements(bdt)} )
	print('entanglements')
	print(etime)
	e <- heights_entanglements[[2]]
	
	pltime <- system.time( { ll <- calculate.coalescent.log.pseudolikelihood.0(sampleTime, bdt, heights_entanglements, FGY=FGY, integrationMethod = 'rk4', discretizeRates=TRUE, fgyResolution = 200) } )
	print(ll)
	print(pltime)
}

sample.initial.conditions <- function()
{
	theta <- c( 
	  runif(3, .05, .5)
	  , runif(1, .5, 2) )
	names(theta) <- ESTNAMES
	as.list( log(theta) )
}

#~ fit the model
if (TRUE)
{
#~ 	X11(); plot(times0, Ys[, 'I0'], type='lines', ylim = c(0, 800), main='true trajectories')
#~ 		for (n in INFECTEDNAMES){ lines(times0, Ys[,n], lty='solid' ) }
	
	.obj2 <- function(thet)
	{ # base coalescent likelihood
		theta1 <- parms_truth
		theta1[ESTNAMES] <- exp(thet)
		print(unlist(theta1[ESTNAMES]))
		FGY <- solve.ode.model.set.fgy(theta1)
		-coalescent.log.likelihood(bdt, FGY=FGY, integrationMethod = 'euler', finiteSizeCorrections=FALSE, maxHeight=20, discretizeRates=TRUE, fgyResolution = 200)
	}
	
	etime <- system.time( {heights_entanglements <- calculate.coalescent.entanglements(bdt)} )
	print('entanglements')
	print(etime)
	
	.obj3 <- function(beta0, beta1, beta2, I0_0)
	{ #pseudolikelihood (entanglements)
		theta1 <- parms_truth
		theta1[ESTNAMES] <- exp(c(beta0, beta1, beta2, I0_0))
		print(unname( unlist(theta1[ESTNAMES])) )
		FGY <- solve.ode.model.set.fgy(theta1)
		mll <- -calculate.coalescent.log.pseudolikelihood.0(sampleTime, bdt, heights_entanglements, FGY=FGY, integrationMethod = 'rk4', discretizeRates=TRUE, fgyResolution = 200)
		#round(mll, 3)
		mll
	}
	
	print(paste(date(), 'fitting model'))
	
	NFITS <- 10
	of <- .obj3 #.obj4 #.obj3
	
	
	
	
	#FITS 
	fits <- lapply( 1:NFITS, function(i)
	{
		print(paste(date(), 'start fit', i))
		of0 <- Inf
		while(is.infinite(of0) || is.na(of0)) 
		{
			theta <- sample.initial.conditions()
			of0 <- do.call( of, theta )
		}
		print(paste('initial conditions:', exp(unlist(theta))))
		#tryCatch( { optim(theta, of,  control=list(trace=6,reltol=1e-04,maxit=300)) }, 
		#  error = function(e) browser() )
		fit <- mle2(of, start = as.list(theta), method='Nelder-Mead', skip.hessian=TRUE, control=list(trace=6,reltol=1e-04,maxit=300,  hessian=FALSE) )
	})
	print(paste(date(), 'done with fits iter 1'))
	bestfit_i <- which.min( sapply(1:NFITS, function(i) -logLik(fits[[i]])) )
	fit <- fits[[bestfit_i]]
	if (attr(fit, 'details')$convergence!=0)
	{
		fit = fit_ml <- mle2(of, start = as.list(coef(fit)), method='Nelder-Mead', control=list(trace=6,reltol=1e-04,maxit=300,  hessian=FALSE) )
		# optim(fit$par, of,  control=list(trace=6,reltol=1e-04,maxit=300))
	}
	
	
	
	
	#OUTPUT
	print( summary(fit) )
	results <- rbind( unlist(parms_truth[ESTNAMES])
	  , exp(coef(fit))
	 )
	colnames(results) <- ESTNAMES
	print(results)
	
	#~  print all fits
	print('')
	allfits <- t(
	  sapply(fits, function(f) c( logLik(f) , attr(f,'details')$convergence, exp( coef(f) ) ) )
	)
	colnames(allfits) <- c('val', 'convergence', ESTNAMES)
	af_i <- sort( allfits[,'val'], index.return=TRUE)$ix
	allfits <- allfits[af_i,]
	print ( allfits )
	
	# approximate confidence intervals
	print('')
	cis <- unname( rbind( exp( coef(fit) - diag( sqrt( vcov(fit) )) * 1.96 ), 
	 exp( coef(fit) + diag( sqrt( vcov(fit) )) * 1.96 ) ) )
	 colnames(cis) <- ESTNAMES
	print( cis )
	
	
	
	print('')
	print('obj function at true_parms')
	print(do.call(of, as.list(theta0)) )
	print('obj function at fitted_parms')
	print(do.call(of, as.list(coef(fit))) )
	
	
	
	#PLOTS 
	# plot MPL values
	X11(); plot( sort(sapply(fits, function(f) logLik(f) )  ) , main='MPL values for all fits')
	#sapply(fits, function(f) c( exp( f$par ), f$value , f$convergence) )
	
	# plot num infected
		X11(); plot(times0, Ys[, 'I0'], type='lines', ylim = c(0, 800), main='MPLE')
		for (n in INFECTEDNAMES){ lines(times0, Ys[,n], lty='solid' ) }
		parameters <- parms_truth
		parameters[ESTNAMES] <- exp(coef(fit))
			FGY <- solve.ode.model.set.fgy(parameters)
			o <- t(sapply(times0, FGY$Y.))
			colnames(o) <- INFECTEDNAMES
			for (n in INFECTEDNAMES) { lines( times0,  o[,n], lty='dashed' ) }
	#
	
	# plot prop transmissions from each stage
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
		plot(heights, true_paf[, 'I0'], type='lines', ylim = c(0, 1), main='MPLE')
		for (n in INFECTEDNAMES){ lines(heights, true_paf[,n], lty='solid', col=pcol[n] ) }
		fitted_paf <- paf.through.time( coef(fit) )
		for (n in INFECTEDNAMES) { lines( heights,  fitted_paf[,n], lty='dashed' , col=pcol[n]) }
		print(true_paf[1,])
		print(fitted_paf[1,])
	
	
		
}


if (FALSE)
{
	#  likelihood profiles
	print(paste(date(), 'start profile'))
	# seems like std.err need to be manually supplied for this to work
	proftime <- system.time( { prof_beta0 <- profile(fit, which=1,tol.newmin=Inf, maxsteps=50, trace=TRUE, std.err=.1) })
	print(paste(date(), 'end profile'))
	print(proftime)
	print( exp( confint(prof_beta0) ))
	X11()
	plot(prof_beta0)
	# 1.5 minutes to do one parameter... but maybe not reliable
	#profile(fit, which = 1, maxsteps = 10, alpha = 0.05,trace =TRUE)
	
	proftime <- system.time( { prof_beta1 <- profile(fit, which=2,tol.newmin=Inf, maxsteps=50, trace=TRUE, std.err=.1) })
	print(paste(date(), 'end profile'))
	print(proftime)
	print( exp( confint(prof_beta1) ))
	X11()
	plot(prof_beta1)
	
	proftime <- system.time( { prof_beta2 <- profile(fit, which=3,tol.newmin=Inf, maxsteps=50, trace=TRUE, std.err=.1) })
	print(paste(date(), 'end profile'))
	print(proftime)
	print( exp( confint(prof_beta2) ))
	X11()
	plot(prof_beta2)
}

#~ 1.7 sec to calc stats from tree
#~ .5 sec to simulate tree
#~ 2.2 total
#~ .5 sec to analytically solve moments for model
