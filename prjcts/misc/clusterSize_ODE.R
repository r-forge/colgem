################################################
#
#	script helper functions
#
################################################

# comparison of model and empirical trees
colgem.project.moments.for.odesystem.comparisonplots <- function(eM, mM)
{
	X11()
	par(mfrow=c(3,3))
	for (i in 1:3)
	{
		for (j in 1:3)
		{
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


colgem.project.moments.for.odesystem<- function()
{
	require(deSolve)
	
	#~ parms_truth <<- list(gamma0 = 1, gamma1 = 1/7, gamma2 = 1/2, mu = 1/30, b=.036, beta0 = 0.775, beta1=0.08, beta2=0.08, S0=2500, alpha = .05) 
	parms_truth 	<<- list(gamma0 = 1, gamma1 = 1/7, gamma2 = 1/2, mu = 1/30, b=.036, beta0 = 1+1/30, beta1=(1+1/30)/10, beta2=(1+1/30)/10, S0=5000, alpha = 4) 
	FGYPARMS 		<<- parms_truth
	INFECTEDNAMES 	<<- c('I0', 'I1', 'I2')
	
	################################################
	#
	#	define F G solve.model.set.fgy that are assumed in rcolgem.R
	#
	################################################	
	#	infection of S into I0
	F. <<- function(t, x = NA)
	{
		if (is.na(x[1]) )  
			x 	<- x.interp(t) 
		N 		<- sum(x) 
		I 		<- sum(x[INFECTEDNAMES])
		lambda	<- exp(-FGYPARMS$alpha * I/N) * c( FGYPARMS$beta0*x['I0'], FGYPARMS$beta1*x['I1'] , FGYPARMS$beta2*x['I2'])  * x['S']/N
		unname( cbind( lambda, matrix(0, nrow=3, ncol=2) ) )
	}	
	#	state transitions flows between model compartments
	#	NOT USED
	G. <<- function(t, x = NA)			
	{
		if (is.na(x[1]) )  
			x 	<- x.interp(t) 
		g 		<- matrix(0, nrow=3, ncol=3)
		g[1,2] 	<- FGYPARMS$gamma0 * x['I0']
		g[2,3] 	<- FGYPARMS$gamma1 * x['I1']
		g
	}	
	#	define ODE 
	dx <- function(t, x, parms, ...) 
	{
		with(parms, {
					lambda 	<- sum( F.(t, x = x))
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
	solve.model.set.fgy <- function(parameters)
	{ 
		FGYPARMS 	<<- parameters
		times 		<- seq(from=0, to=50, length.out = 100)
		x0 			<- c(S = FGYPARMS$S0, I0=1, I1=0, I2=0) 
		o 			<- ode( x0, times, func=dx, parms=parameters, method='adams')
		
		x.interps	<-  sapply( names(x0), function(n) approxfun(o[,1], o[,n], yleft=o[1,n], yright=o[nrow(o),n]) )
		x.interp 	<<- function(t) { sapply(x.interps, function(interp) interp(t) ) }  
		Y. 			<<- function(t) x.interp(t)[INFECTEDNAMES]
	}	
	
	################################################
	#
	#	start script
	#
	################################################
	
	# simulate coalescent tree with true parameters 
	solve.model.set.fgy(parms_truth) 
	phi				<- .50 # sample fraction
	sampleTime 		<- 50
	Y.sampleTime 	<- Y.(sampleTime)
	stateIndices 	<- rep( 1:3, round( phi * Y.sampleTime ) ) # sample each of three states in proportion to size
	n 				<- length(stateIndices)
	sampleTimes 	<- rep(sampleTime, n)
	sampleStates 	<- diag(3)[stateIndices,]
	bdt 			<- simulatedBinaryDatedTree(sampleTime, sampleStates, discretizeRates=TRUE) 
	
	# calculate empirical stats
	heights 		<- seq(1, 50, length.out=50)
	empiricalMoments<- eM	<- calculate.cluster.size.moments.from.tree(bdt, heights)
	
	# calculate model stats
	modelMoments	<- mM 	<- calculate.cluster.size.moments.from.model(sampleTime, sampleStates , timeResolution = 50, discretizeRates=TRUE, fgyResolution = 100 , integrationMethod = 'rk4')
	colgem.project.moments.for.odesystem.comparisonplots(eM, mM) 
	
	#~ now do a comparison with model parameters different from true parameters 
	parameters 			<- parms_truth
	parameters$beta0 	<- parms_truth$beta0/2
	parameters$beta1 	<- parms_truth$beta1*5
	solve.model.set.fgy(parameters)
	modelMoments 		<- mM 	<- calculate.cluster.size.moments.from.model(sampleTime, sampleStates , timeResolution = 50, discretizeRates=TRUE, fgyResolution = 100 , integrationMethod = 'rk4')
	colgem.project.moments.for.odesystem.comparisonplots(eM, mM) 
}

