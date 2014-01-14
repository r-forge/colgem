################################################
#
#	script helper functions
#
################################################

# comparison of model and empirical trees
cg.sde.comparison.plots <- function(eM, mM)
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

cg.sde.comparison.plots2 <- function(eM, nsims, solve.model.set.fgy, parms)
{
	mM <- list()
	for (isim in 1:nsims)
	{
		mM_time<- system.time({  
					modelMoments = mM_i <- calculate.cluster.size.moments.from.model(sampleTime, sampleStates , timeResolution = 50, discretizeRates=TRUE, fgyResolution = 100 , integrationMethod = 'rk4')
				})
	}
	mM <- lapply(1:nsims, function(isim) 
			{
				solve.model.set.fgy(parms) ;
				calculate.cluster.size.moments.from.model(sampleTime, sampleStates , timeResolution = 50, discretizeRates=TRUE, fgyResolution = 100 , integrationMethod = 'rk4');
			})
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

cg.sde.modelMoments<- function(parms, eM, sampleTime, sampleStates, solve.model.set.fgy, nsims)
{
	mM <- list()
	for (isim in 1:nsims)
	{
		mM.time <- system.time({  
					modelMoments = mM_i <- calculate.cluster.size.moments.from.model(sampleTime, sampleStates , timeResolution = 50, discretizeRates=TRUE, fgyResolution = 100 , integrationMethod = 'rk4')
				})
	}
	mM <- lapply(1:nsims, function(isim) 
			{
				solve.model.set.fgy(parms)
				calculate.cluster.size.moments.from.model(sampleTime, sampleStates , timeResolution = 50, discretizeRates=TRUE, fgyResolution = 100 , integrationMethod = 'rk4');
			})
	list(mM=mM, mM.time=mM.time)	
}

cg.sde.comparison.plots3<- function(heights, eM, mM, dir.name)
{ 
	nsims	<- length(mM)
	# plot deltas
	file	<- paste(dir.name,'delta_Mij.pdf', sep='/')
	cat(paste("plot deltas to file=",file))
	pdf( file )
	par( mfrow=c(3,3) )
	nheights 	<- length(heights)
	for (i in 1:3)
	{
		for (j in 1:3)
		{
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
	
	# plot residuals	
	file	<- paste(dir.name,'residuals_Mij.pdf', sep='/')
	cat(paste("plot residuals to file=",file))
	pdf( file )
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
	file	<- paste(dir.name,'residualshist_Mij.pdf', sep='/')
	cat(paste("plot residuals to file=",file))
	pdf( file )
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

#'	define F G solve.model.set.fgy that are assumed in rcolgem.R
cg.sde.define<- function()
{
	#~ the 'skeleton' functions are deterministic functions of system state
	F.skeleton<- function(t, x = NA)
	{
		if (is.na(x[1]) )  x <- x.interp(t) 
		N <- sum(x) 
		I <- sum(x[INFECTEDNAMES])
		lambda <- exp(-FGYPARMS$alpha * I/N) * c( FGYPARMS$beta0*x['I0'], FGYPARMS$beta1*x['I1'] , FGYPARMS$beta2*x['I2'])  * x['S']/N
		unname( cbind( lambda, matrix(0, nrow=3, ncol=2) ) )
	}
	G<- function(t, x = NA)
	{
		if (is.na(x[1]) )  x <- x.interp(t) 
		g <- matrix(0, nrow=3, ncol=3)
		g[1,2] <- FGYPARMS$gamma0 * x['I0']
		g[2,3] <- FGYPARMS$gamma1 * x['I1']
		g
	}	
	#~ solve using discrete euler method	
	step.x<- function(t, x, parms ) 
	{
		timestep	<- diff(parms$times)[1]
		with(parms, {
					lambda 				<- F.skeleton(t, x = x)
					lambdaTimestep 		<- lambda * timestep
					sumlambda 			<- sum(lambda) 
					sumlambdaTimestep 	<- sum(lambdaTimestep) 
					rsumlambda 			<- abs( sumlambda*timestep + sqrt(sumlambda) * rnorm(1, mean=0, sd=sqrt(timestep)) )
					lambda 				<- lambdaTimestep * rsumlambda/sumlambdaTimestep
					N 					<- sum(x)
					list( x = x + c(
									S = -rsumlambda + timestep*b * N - timestep * mu * x['S'], 
									I0 = rsumlambda - timestep * (mu+gamma0)*x['I0'],
									I1 = timestep*gamma0*x['I0'] - timestep*(mu+gamma1)*x['I1'],
									I2 = timestep*gamma1*x['I1'] - timestep*(mu+gamma2)*x['I2']), 
							lambda = lambda )
				})
	}
	solve.model.set.fgy<- function(parameters, file=NA)
	{ 
		FGYPARMS	<<- parameters		
		lambdas 	<- list()
		x 			<-  c(S = parameters$S_0, I0=parameters$I0_0, I1=parameters$I1_0, I2=parameters$I2_0)
		times		<- parameters$times
		timestep	<- diff(parameters$times)[1]
		m			<- 3
		o 			<- c(0, x)
		for (itimes in 2:length(times))
		{
			t 					<- times[itimes]
			x_lambda 			<- step.x( t, x, parameters)
			x 					<- x_lambda[['x']]
			lambdas[[itimes]] 	<- x_lambda[['lambda']]
			o					<- rbind(o, c(t, x))
		}
		x.interps	<-  sapply( names(x), function(n) approxfun(o[,1], o[,n], yleft=o[1,n], yright=o[nrow(o),n]) )
		x.interp 	<<- function(t) { sapply(x.interps, function(interp) interp(t) ) }  
		Y. 			<<- function(t) x.interp(t)[INFECTEDNAMES]
		F.interps 	<-  list()
		for (k in 1:parameters$m)
		{ 
			F.interps[[k]] 			<- list()
			for (l in 1:parameters$m)
			{
				lambdaskl 			<- sapply(2:length(times), function(itime) lambdas[[itime]][k,l])
				F.interps[[k]][[l]] <- approxfun( times[1:(length(times)-1)], lambdaskl/timestep, method='constant', rule=2)
			}
		}
		F. 			<<- function(t.) t( sapply(1:parameters$m, function(k)  sapply(1:parameters$m, function(l) F.interps[[k]][[l]](t.)) ) )
		if(!is.na(file))
		{
			class(o) <- 'deSolve'
			cat(paste("\nplot to file=",file))
			pdf( file= file, 5, 5 )
			plot(o)
			dev.off()
		}
	}
	
	list(F.skeleton=F.skeleton, G=G, step.x=step.x, solve.model.set.fgy=solve.model.set.fgy)
}

cg.sde<- function()
{
	#	for a single model simulation, compute nsim=20 moment trajectories. 
	#	the randomness comes from fluctuations in the state variables S I0 I1 I2 
	#cg.sde.nsim.mM()
	
	cg.sde.varyparam()
}

cg.sde.varyparam<- function()
{
	my.mkdir(HOME, 'MOMSDE' )
	dir.name		<- paste(HOME, 'MOMSDE',sep='/')
	#define model
	tmp					<- cg.sde.define()
	F.skeleton			<<- tmp$F.skeleton
	G.					<<- tmp$G
	step.x				<<- tmp$step.x
	solve.model.set.fgy	<- tmp$solve.model.set.fgy
	
	#parameter template
	parms.template 	<- list(	m=3, gamma0 = 1, gamma1 = 1/7, gamma2 = 1/2, mu = 1/30, b=.036, beta0 = 1+1/30, beta1=(1+1/30)/10, beta2=(1+1/30)/10, 
								S_0=5000, I0_0=1, I1_0=1, I2_0=1, alpha = 4,
								times=seq(0, 50, by=.1))	
	INFECTEDNAMES 	<<- c('I0', 'I1', 'I2')
	phi				<- .50 # sample fraction
	sampleTime 		<- 50
	nsims			<- 20
							
	#bits of the model that are varied
	gi			<- c(0.25, 0.5, 1, 3) 		#avg duration of I0
	bf			<- c(2,4,8,16) 				#dampening of beta0
	parms.vary	<- expand.grid(gi=gi, bf= bf)	
	parms.vary	<- cbind(parms.vary, gamma0= 1/parms.vary[,'gi'], beta1= parms.template[['beta0']]/parms.vary[,'bf'])
	
	#plot resulting epidemics 
	dummy		<- lapply(seq_len(nrow(parms.vary)),function(i)
			{					  		
				parms.template[c('gamma0','beta1','beta2')]	<- parms.vary[i,c('gamma0','beta1','beta1')]
				
				parms_truth		<<- parms.template	
				FGYPARMS 		<<- parms_truth
				file			<- paste("varyparam.gi=",parms.vary[i,'gi'],"_bf=",parms.vary[i,'bf'],".pdf",sep='')
				file			<- paste(dir.name,file,sep='/')
				solve.model.set.fgy(parms_truth, file)				
			})	
	#with increasing bf: 	smaller depletion of Susc, epidemic more stationary for bf=16  
	
}

cg.sde.nsim.mM<- function()
{	
	my.mkdir(HOME, 'MOMSDE' )
	dir.name		<- paste(HOME, 'MOMSDE',sep='/')
	
	#~ parms_truth <<- list(gamma0 = 1, gamma1 = 1/7, gamma2 = 1/2, mu = 1/30, b=.036, beta0 = 0.775, beta1=0.08, beta2=0.08, S0=2500, alpha = .05) 
	parms_truth 	<<- list(	m=3, gamma0 = 1, gamma1 = 1/7, gamma2 = 1/2, mu = 1/30, b=.036, beta0 = 1+1/30, beta1=(1+1/30)/10, beta2=(1+1/30)/10, 
								S_0=5000, I0_0=1, I1_0=1, I2_0=1, alpha = 4,
								times=seq(0, 50, by=.1)) 
	FGYPARMS 		<<- parms_truth
	INFECTEDNAMES 	<<- c('I0', 'I1', 'I2')
	phi				<- .50 # sample fraction
	sampleTime 		<- 50
	nsims			<- 20
	
	tmp					<- cg.sde.define()
	F.skeleton			<<- tmp$F.skeleton
	G.					<<- tmp$G
	step.x				<<- tmp$step.x
	solve.model.set.fgy	<- tmp$solve.model.set.fgy
	
	################################################
	#
	#	start script
	#
	################################################
	
	
	# simulate coalescent tree with true parameters
	file			<- paste(dir.name,'nsim.mM.epidemic.pdf',sep='/')
	solve.model.set.fgy(parms_truth, file) 
	Y.sampleTime 	<- Y.(sampleTime)
	m				<- 3
	stateIndices 	<- rep( 1:m, round( phi * Y.sampleTime ) ) # sample each of three states in proportion to size
	n 				<- length(stateIndices)
	sampleTimes 	<- rep(sampleTime, n)
	sampleStates 	<- diag(m)[stateIndices,]
	cat(paste("\nsimulate structured coalescent tree, including tip states (0/1)"))
	coalescentTree_time <- system.time({
				bdt <- simulatedBinaryDatedTree(sampleTime, sampleStates, discretizeRates=TRUE) 
			})
	# calculate empirical stats
	heights 		<- seq(0, 50, length.out=50)
	cat(paste("\ncalculate moments of observed tree by averaging over lineages; tip states 0/1"))
	empiricalMoments= eM <- calculate.cluster.size.moments.from.tree(bdt, heights)
	cat(paste("\ncalculate moments under model"))
	tmp				<- cg.sde.modelMoments(parms_truth, eM, sampleTime, sampleStates, solve.model.set.fgy, nsims)
	mM				<- tmp$mM
	cat(paste("\nplot results"))
	cg.sde.comparison.plots3(heights, eM , mM, dir.name)	
}


# calculate model stats
#~ mM_time <- system.time({
#~   modelMoments = mM <- calculate.cluster.size.moments.from.model(sampleTime, sampleStates , timeResolution = 50, discretizeRates=TRUE, fgyResolution = 100 , integrationMethod = 'rk4')
#~ })
#~ cg.sde.comparison.plots(eM, mM) 

#~ now do a comparison with model parameters different from true parameters 
#~ parameters <- parms_truth
#~ parameters$beta0 <- parms_truth$beta0/2
#~ parameters$beta1 <- parms_truth$beta1*5
#~ solve.model.set.fgy(parameters)
#~ modelMoments = mM <- calculate.cluster.size.moments.from.model(sampleTime, sampleStates , timeResolution = 50, discretizeRates=TRUE, fgyResolution = 100 , integrationMethod = 'rk4')
#~ cg.sde.comparison.plots(eM, mM, solve.model.set.fgy, parms_truth) 


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
