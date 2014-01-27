################################################
#
#	script helper functions
#
################################################################################################

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
################################################################################################
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
################################################################################################
cg.sde.modelMoments<- function(parms, sampleTime, sampleStates, solve.model.set.fgy, nsims)
{
	ans <- lapply(1:nsims, function(isim) 
			{
				solve.model.set.fgy(parms)
				#timeResolution = 50; discretizeRates=TRUE; fgyResolution = 100; integrationMethod = 'rk4'
				calculate.cluster.size.moments.from.model(sampleTime, sampleStates , timeResolution = 50, discretizeRates=TRUE, fgyResolution = 100 , integrationMethod = 'rk4');
			})
	ans
}
################################################################################################
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
################################################################################################
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
		list(simu=o, x.interp=x.interp, Y=Y., F=F.)
	}
	
	list(F.skeleton=F.skeleton, G=G, step.x=step.x, solve.model.set.fgy=solve.model.set.fgy)
}
################################################################################################
util.fade.col<-function(col,alpha=0.5)
{
	return(rgb(col2rgb(col)[1]/255,col2rgb(col)[2]/255,col2rgb(col)[3]/255,alpha))
}
################################################################################################
cg.sde<- function()
{	
	require(data.table)
	#	for a single model simulation, compute nsim=20 moment trajectories. 
	#	the randomness comes from fluctuations in the state variables S I0 I1 I2 
	#cg.sde.nsim.mM()
	
	#	run one simulation of the epidemic under the cg.sde model for various parameters
	#cg.sde.varyparam()
	
	#	produce 1e2 pseudo data sets under the cg.sde model for various paramaters
	#cg.sde.get.pseudodata()
	
	#	evaluate different pseudo likelihoods to the empirical cluster size distributions
	#cg.sde.eval.pseudolkl.to.empirical.clustersize.distribution()
	
	#	evaluate different pseudo likelihoods to the eMis
	cg.sde.eM.pseudolkl()
}
################################################################################################
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
								times=seq(0, 50, by=.1), sampleTime=50, phi=0.5)	
	INFECTEDNAMES 	<<- c('I0', 'I1', 'I2')
	nsims			<- 20
							
	#bits of the model that are varied
	gi			<- c(0.5, 1, 3) 			#avg duration of I0
	bf			<- c(2,4,8,16) 				#dampening of beta0
	parms.vary	<- expand.grid(gi=gi, bf= bf)	
	parms.vary	<- cbind(parms.vary, gamma0= 1/parms.vary[,'gi'], beta1= parms.template[['beta0']]/parms.vary[,'bf'])
	
	#plot resulting epidemics for varying parameters
	dummy		<- lapply(seq_len(nrow(parms.vary)),function(i)
			{					  		
				parms.template[c('gamma0','beta1','beta2')]	<- parms.vary[i,c('gamma0','beta1','beta1')]
				
				parms_truth		<<- parms.template	
				FGYPARMS 		<<- parms_truth
				file			<- paste("varyparam.gi=",parms.vary[i,'gi'],"_bf=",parms.vary[i,'bf'],".pdf",sep='')
				file			<- paste(dir.name,file,sep='/')
				ans				<- solve.model.set.fgy(parms_truth, file)				
			})	
	#with increasing bf: 	smaller depletion of Susc, epidemic more stationary for bf=16  
	#with increasing gi:	depletion of susceptibles with gi~3. epidemic does not take off with gi~0.25, need at least 0.5
}
################################################################################################
cg.sde.get.pseudodata.for.param<- function(parms, bdt.n= 1e2, bdt.i=NA, bdt.heights=seq(0, 50, length.out=50), distr.heights=seq(10, 50, length.out=50), file=NA)
{
	#define model
	INFECTEDNAMES 		<<- c('I0', 'I1', 'I2')
	tmp					<-  cg.sde.define()
	F.skeleton			<<- tmp$F.skeleton
	G.					<<- tmp$G
	step.x				<<- tmp$step.x
	solve.model.set.fgy	<-  tmp$solve.model.set.fgy
	
	#set up parameters
	m				<-  parms$m
	parms_truth		<<- parms	
	FGYPARMS 		<<- parms_truth
	
	if(!is.na(file))
	{
		if( substr(file, nchar(file)-1,nchar(file))!='.R') stop("file expected to end in .R")
		file	<- substr(file, 1, nchar(file)-2)
	}
	if(is.na(bdt.i))
		bdt.i	<- seq_len(bdt.n)
	#simulated btd.n pseudo data sets
	pseudo.datasets	<- lapply(bdt.i, function(i)
			{
				if(!is.na(file))
				{
					file.bdt					<- paste(file,'_bdtn',i,'.R',sep='')
					cat(paste("\ntry load ",file.bdt))
					options(show.error.messages = FALSE, warn=1)		
					readAttempt					<-try(suppressWarnings(load(file.bdt)))						
					options(show.error.messages = TRUE)								
				}	
				if(inherits(readAttempt, "try-error"))
				{	
					cat(paste('\nprocess bdt.i',i))
					pseudo.data					<- solve.model.set.fgy(parms_truth)	
					pseudo.data$stateIndices 	<- rep( 1:m, round( parms$phi * pseudo.data$Y( parms$sampleTime ) ) ) # sample each of three states in proportion to size	 
					pseudo.data$sampleTimes 	<- rep(parms$sampleTime, length(pseudo.data$stateIndices) )
					pseudo.data$sampleStates 	<- diag(m)[pseudo.data$stateIndices,]
					pseudo.data$bdt 			<- simulatedBinaryDatedTree(parms$sampleTime, pseudo.data$sampleStates, discretizeRates=TRUE)				
					pseudo.data$distr			<- calculate.cluster.size.distr.from.tree(pseudo.data$bdt, distr.heights)
					tmp			 				<- calculate.cluster.size.moments.from.tree(pseudo.data$bdt, bdt.heights)				
					pseudo.data$eMi 			<- tmp$Mi
					pseudo.data$eMij 			<- tmp$Mij				
					tmp							<- lapply( seq_len(parms$mM.replicate), function(b)
						{
							#	get moments for new model solution to F G Y
							solve.model.set.fgy(parms)				
							tmp					<- calculate.cluster.size.moments.from.model(parms$sampleTime, pseudo.data$sampleStates , timeResolution = 50, discretizeRates=TRUE, fgyResolution = 100 , integrationMethod = 'rk4')
							#	store moments as data.table
							tmp$Mi				<- cbind( as.data.table(t( tmp$Mi )), height=tmp$heights, replicate=b )				
							tmp2				<- sapply( seq_len( dim(tmp$Mij)[3] ), function(h)	tmp$Mij[,,h][upper.tri(tmp$Mij[,,1],diag=1)] )
							rownames(tmp2)		<- sapply(1:m, function(i) paste('state.',i,1:m,sep=''))[upper.tri(tmp$Mij[,,1],diag=1)] 
							tmp$Mij				<- cbind( as.data.table(t( tmp2 )), height=tmp$heights, replicate=b )
							tmp
						})
					pseudo.data$mMi				<- lapply( seq_along(tmp), function(b)	tmp[[b]]$Mi	)
					pseudo.data$mMi				<- do.call('rbind',pseudo.data$mMi	)
					pseudo.data$mMij			<- lapply( seq_along(tmp), function(b)	tmp[[b]]$Mij	)
					pseudo.data$mMij			<- do.call('rbind',pseudo.data$mMij	)
					if(!is.na(file))
					{					
						cat(paste('\nsave pseudo.data to',file.bdt))
						save(pseudo.data, file=file.bdt)
					}	
				}								
				pseudo.data		
			})
	pseudo.datasets	
}
################################################################################################
cg.sde.get.mM.for.param<- function(parms, file=NA)
{
	require(data.table)
	#define model
	INFECTEDNAMES 		<<- c('I0', 'I1', 'I2')
	tmp					<-  cg.sde.define()
	F.skeleton			<<- tmp$F.skeleton
	G.					<<- tmp$G
	step.x				<<- tmp$step.x
	solve.model.set.fgy	<-  tmp$solve.model.set.fgy
	
	#set up parameters
	m				<-  parms$m
	parms_truth		<<- parms	
	FGYPARMS 		<<- parms_truth
	
	if(!is.na(file))
	{
		if( substr(file, nchar(file)-1,nchar(file))!='.R') stop("file expected to end in .R")
		cat(paste("\ntry load ",file))
		options(show.error.messages = FALSE, warn=1)		
		readAttempt					<-try(suppressWarnings(load(file)))						
		options(show.error.messages = TRUE)								
	}
	if(is.na(file) | inherits(readAttempt, "try-error"))
	{		
		tmp			<- lapply( seq_len(parms$mM.replicate), function(b)
				{
					#	get moments for new model solution to F G Y
					pseudo.data					<- solve.model.set.fgy(parms_truth)	
					pseudo.data$stateIndices 	<- rep( 1:m, round( parms$phi * pseudo.data$Y( parms$sampleTime ) ) ) # sample each of three states in proportion to size	 
					pseudo.data$sampleTimes 	<- rep(parms$sampleTime, length(pseudo.data$stateIndices) )
					pseudo.data$sampleStates 	<- diag(m)[pseudo.data$stateIndices,]					
					#
					tmp					<- calculate.cluster.size.moments.from.model(parms$sampleTime, pseudo.data$sampleStates , timeResolution = 50, discretizeRates=TRUE, fgyResolution = 100 , integrationMethod = 'rk4')
					#	store moments as data.table
					tmp$Mi				<- cbind( as.data.table(t( tmp$Mi )), height=tmp$heights, replicate=b )				
					tmp2				<- sapply( seq_len( dim(tmp$Mij)[3] ), function(h)	tmp$Mij[,,h][upper.tri(tmp$Mij[,,1],diag=1)] )
					rownames(tmp2)		<- sapply(1:m, function(i) paste('state.',i,1:m,sep=''))[upper.tri(tmp$Mij[,,1],diag=1)] 
					tmp$Mij				<- cbind( as.data.table(t( tmp2 )), height=tmp$heights, replicate=b )
					tmp
				})
		df.mMi		<- lapply( seq_along(tmp), function(b)	tmp[[b]]$Mi	)
		df.mMi		<- do.call('rbind',df.mMi	)
		df.mMi.raw	<- copy( df.mMi )
		df.mMij		<- lapply( seq_along(tmp), function(b)	tmp[[b]]$Mij	)
		df.mMij		<- do.call('rbind',df.mMij	)				
		
		#
		#	clean up
		#
		
		#	clean up df.mMi: remove numerical inaccuracies in mMi
		cmd		<- paste('is.nan(state.',1:m,')',' | ','abs(state.',1:m,')<',1e-6,sep='', collapse=' | ')	
		cmd		<- paste('df.mMi[, which(',cmd,')]',sep='')
		tmp		<- eval( parse(text=cmd) )
		set(df.mMi, tmp, 'height', NA)				
		cmd		<- paste('{ ok.last<- which(is.na(height));  ok.last<- ifelse(length(ok.last), ok.last-1, length(height));    list(', paste('state.',1:m,'= state.',1:m,'[seq_len(ok.last)]',sep='', collapse=','),', height=height[seq_len(ok.last)]) }', sep='')
		cmd		<- paste('df.mMi[,',cmd,', by="replicate"]')
		df.mMi	<- eval( parse(text=cmd) )
		
		#	remove negative changes in mMi -- ONLY FOR HOMOCHRONEOUS SAMPLING AND NONDECREASING POPSIZE
		cmd		<- paste('state.',1:m,'<=',0,sep='', collapse=' | ')	
		cmd		<- paste('df.mMi[, which(',cmd,')]',sep='')
		tmp		<- eval( parse(text=cmd) )
		set(df.mMi, tmp, 'height', NA)				
		cmd		<- paste('{ ok.last<- which(is.na(height));  ok.last<- ifelse(length(ok.last), ok.last-1, length(height));    list(', paste('state.',1:m,'= state.',1:m,'[seq_len(ok.last)]',sep='', collapse=','),', height=height[seq_len(ok.last)]) }', sep='')
		cmd		<- paste('df.mMi[,',cmd,', by="replicate"]')
		df.mMi	<- eval( parse(text=cmd) )
		
		#
		#	get delta mMi
		#
				
		df.mMi[, dummy:= -height]		
		setkey(df.mMi, replicate, dummy)
		cmd		<- paste('state.',1:m,'= -diff(state.',1:m,')',sep='', collapse=',')
		cmd		<- paste('df.mMi[, list(',cmd,', height=height[-1]), by="replicate"]',sep='')
		df.mMid	<- eval( parse(text=cmd) )
		setkey(df.mMid, replicate, height)
		
		#	clean up df.mMid: remove numerical inaccuracies 
		cmd		<- paste('state.',1:m,'<',0,sep='', collapse=' | ')	
		cmd		<- paste('df.mMid[, which(',cmd,')]',sep='')
		tmp		<- eval( parse(text=cmd) )
		set(df.mMid, tmp, 'height', NA)	
		cmd		<- paste('{ ok.last<- which(is.na(height));  ok.last<- ifelse(length(ok.last), ok.last-1, length(height));  list(', paste('state.',1:m,'= state.',1:m,'[seq_len(ok.last)]',sep='', collapse=','),', height=height[seq_len(ok.last)]) }', sep='')
		cmd		<- paste('df.mMid[,',cmd,', by="replicate"]')
		df.mMid	<- eval( parse(text=cmd) )
		
		
		sim			<- list()
		sim$parms	<- parms			
		sim$mMi		<- df.mMi
		sim$mMi.raw	<- df.mMi.raw
		sim$mMid	<- df.mMid
		sim$mMij	<- df.mMij										
	}
	if(!is.na(file))
	{					
		cat(paste('\nsave sim to',file))
		save(sim, file=file)
	}
	sim	
}
################################################################################################
# plot distribution of cluster sizes for each tree height, 
# this pools over many repeat runs to get a nice distribution
cg.sde.tiptypenumber.pseudolkl.plotdistr<- function(pseudo.datasets, distr, heights, file)
{
	m		<- length( which( grepl('state',colnames(distr)) ) )
	pdf(file=file, w=15, h=5)
	par(mfcol=c(1,m))
	lapply( seq_along(heights), function(ih)
			{
				cat(paste("\nplot distr",ih))
				pseudo_h	<- subset( distr, height==heights[ih] )				
				
				eMi			<- sapply(seq_along(pseudo.datasets), function(k)	pseudo.datasets[[k]][['eMi']][,ih]	)
				eMii		<- sapply(seq_along(pseudo.datasets), function(k)	diag(pseudo.datasets[[k]][['eMij']][,,ih])	)
				eMi.avg		<- rowMeans( eMi )
								
				sapply( seq_len(m), function(j)
						{
							tmp	<- unclass( pseudo_h[,j,with=0] )[[1]]
							hist( tmp, main=paste('state=',j,'height=',round(heights[ih],d=1)), breaks=max(tmp)+1, xlab='counts', col='gray50', bty='n', border=NA )				
							abline(v=eMi.avg[j])	
						})				
			})
	dev.off()
}
################################################################################################
# plot distribution of cluster sizes for each tree height and add possible pseudo likelihoods
# this pools over many repeat runs to get a nice distribution
cg.sde.tiptypenumber.pseudolkl.plotpotentiallkl<- function(pseudo.datasets, distr, file=NA)
{
	library(MASS)
	library(RColorBrewer)
	m			<- length( which( grepl('state',colnames(distr)) ) )
	heights		<- unique(distr[,height])
	plot		<- !is.na(file)
	if(plot)
	{
		pdf(file=file, w=12, h=4)
		par(mfcol=c(1,m))
		cols		<- brewer.pal(3, 'Set1')
		cols		<- rep(cols,each=2)
		ltys		<- c(1,2,1,2,1,2)
		lwds		<- rep(1.5,6)
		legend.txt	<- c('Poisson mle','Poisson mo','Geometric mle','Geometric mo','NegBin mle', 'NegBin mo')
	}
	df.test		<- lapply( seq_along(heights), function(ih)
			{
				cat(paste("\nplot distr",ih))
				pseudo_h	<- subset( distr, height==heights[ih] )				
				
				eMi			<- sapply(seq_along(pseudo.datasets), function(k)	pseudo.datasets[[k]][['eMi']][,ih]	)
				eMii		<- sapply(seq_along(pseudo.datasets), function(k)	diag(pseudo.datasets[[k]][['eMij']][,,ih])	)
				eMi.avg		<- rowMeans( eMi )
				eMii.avg	<- rowMeans( eMii )
				#apply( eMi, 1, sd )
				#apply( eMii, 1, sd )				
				suppressWarnings(
				df.test		<- lapply( seq_len(m), function(j)
						{
							tmp		<- unclass( pseudo_h[,j,with=0] )[[1]]
							tmp		<- unclass( pseudo_h[,j,with=0] )[[1]]				
							hist.x	<- hist( tmp, main=paste('state=',j,'height=',round(heights[ih],d=1)), breaks=max(tmp)+1, xlab='counts', freq=FALSE, col='gray50', bty='n', border=NA, plot=plot )
							hist.y	<- hist.x$density
							fit.x	<- hist.x$breaks[ -length(hist.x$breaks) ]
							
							#	try different pseudo lkls
							#	Poisson fit		
							mle			<- fitdistr(tmp, densfun="Poisson")
							fit.y		<- dpois(fit.x, mle$estimate)
							t.stat		<- sum( ( hist.y - fit.y )^2 / fit.y )
							ans			<- data.table( height=ih, t.stat=t.stat, p= pchisq(t.stat, df=length(fit.x)-2, lower.tail=F), state=j, distr='poi', method='mle' )
							if(plot) lines(fit.x, fit.y, col=cols[1], lty=ltys[1], lwd=lwds[1])
							fit.y		<- dpois(fit.x, eMi.avg[j])
							t.stat		<- sum( ( hist.y - fit.y )^2 / fit.y )
							ans			<- rbind(ans, data.table( height=ih, t.stat=t.stat, p= pchisq(t.stat, df=length(fit.x)-2, lower.tail=F), state=j, distr='poi', method='mo' ))						
							if(plot) lines(fit.x, fit.y, col=cols[2], lty=ltys[2], lwd=lwds[2])
							#	Geometric fit								
							mle			<- fitdistr(tmp, densfun="Geometric")
							fit.y		<- dgeom(fit.x, mle$estimate)
							t.stat		<- sum( ( hist.y - fit.y )^2 / fit.y )
							if(plot) lines(fit.x, fit.y, col=cols[3], lty=ltys[3], lwd=lwds[3])
							ans			<- rbind(ans, data.table( height=ih, t.stat=t.stat, p= pchisq(t.stat, df=length(fit.x)-2, lower.tail=F), state=j, distr='geom', method='mle' ))
							fit.y		<- dgeom(fit.x, 1/(1+eMi.avg[j]))
							t.stat		<- sum( ( hist.y - fit.y )^2 / fit.y )
							ans			<- rbind(ans, data.table( height=ih, t.stat=t.stat, p= pchisq(t.stat, df=length(fit.x)-2, lower.tail=F), state=j, distr='geom', method='mo' ))
							if(plot) lines(fit.x, fit.y, col=cols[4], lty=ltys[4], lwd=lwds[4])
							#	NegBin fit							
							mle			<- fitdistr(tmp, densfun="negative binomial")
							fit.y		<- dnbinom(fit.x, size=mle$estimate['size'], mu=mle$estimate['mu'])
							t.stat		<- sum( ( hist.y - fit.y )^2 / fit.y )
							if(plot) lines(fit.x, fit.y, col=cols[5], lty=ltys[5], lwd=lwds[5])
							ans			<- rbind(ans, data.table( height=ih, t.stat=t.stat, p= pchisq(t.stat, df=length(fit.x)-3, lower.tail=F), state=j, distr='negbin', method='mle' ))
							#	fit Neg Bin based on first two moments
							nbinom.p	<- 1 - eMi.avg[j]/eMii.avg[j]
							nbinom.r	<- eMi.avg[j] * (1 - nbinom.p) / nbinom.p
							fit.y		<- dnbinom(fit.x, size=nbinom.r, prob=1-nbinom.p)
							t.stat		<- sum( ( hist.y - fit.y )^2 / fit.y )
							if(plot) lines(fit.x, fit.y, col=cols[6], lty=ltys[6], lwd=lwds[6])
							ans			<- rbind(ans, data.table( height=ih, t.stat=t.stat, p= pchisq(t.stat, df=length(fit.x)-3, lower.tail=F), state=j, distr='negbin', method='mo' ))																					
							#	plot mean
							if(plot && j==1)	legend('topright', bty='n', border=NA, lty=ltys, fill=cols, legend=legend.txt)
							ans
						})	
				)
				df.test	<- do.call('rbind',df.test)
				df.test				
			})
	if(plot) dev.off()		
	df.test	<- do.call('rbind',df.test)
	df.test			
}
################################################################################################
cg.sde.tiptypenumber.pseudolkl.plotpotentiallkl2<- function(pseudo.datasets, distr, file=NA)
{
	library(MASS)
	library(RColorBrewer)
	m			<- length( which( grepl('state',colnames(distr)) ) )
	heights		<- unique(distr[,height])
	plot		<- !is.na(file)
	if(plot)
	{
		pdf(file=file, w=12, h=4)
		par(mfcol=c(1,m))
		cols		<- brewer.pal(3, 'Set1')
		cols		<- rep(cols,each=2)
		ltys		<- c(1,2,1,2,1,2)
		lwds		<- rep(1.5,6)
		legend.txt	<- c('Poisson mle','Poisson mo','Geometric mle','Geometric mo','NegBin mle', 'NegBin mo')
	}
	df.test		<- lapply( seq_along(heights), function(ih)
			{
				cat(paste("\nprocess height",ih))
				pseudo_h	<- subset( distr, height==heights[ih] )								
				eMi			<- sapply(seq_along(pseudo.datasets), function(k)	pseudo.datasets[[k]][['eMi']][,ih]	)
				eMii		<- sapply(seq_along(pseudo.datasets), function(k)	diag(pseudo.datasets[[k]][['eMij']][,,ih])	)
				#apply( eMi, 1, sd )
				#apply( eMii, 1, sd )				
				suppressWarnings(
						df.test		<- lapply( seq_len(m), function(j)
								{
									tmp		<- unclass( pseudo_h[,j,with=0] )[[1]]
									hist.x	<- hist( tmp, main=paste('state=',j,'height=',round(heights[ih],d=1)), breaks=max(tmp)+1, xlab='counts', freq=FALSE, col='gray50', bty='n', border=NA, plot=plot )
									hist.y	<- hist.x$density
									fit.x	<- hist.x$breaks[ -length(hist.x$breaks) ]									
									#	try different pseudo lkls
									#	Poisson fit		
									mle			<- fitdistr(tmp, densfun="Poisson")
									fit.y		<- dpois(fit.x, mle$estimate)
									if(plot) lines(fit.x, fit.y, col=cols[1], lty=ltys[1], lwd=lwds[1])
									t.stat		<- sum( ( hist.y - fit.y )^2 / fit.y )
									ans			<- data.table( height=ih, t.stat=t.stat, p= pchisq(t.stat, df=length(fit.x)-2, lower.tail=F), state=j, replicate=1, distr='poi', method='mle' )									
									t.stat		<- sapply(seq_along(pseudo.datasets), function(k)
											{
												fit.y		<- dpois(fit.x, eMi[j,k])
												if(plot) lines(fit.x, fit.y, col=util.fade.col(cols[2],0.4), lty=ltys[2], lwd=lwds[2])
												sum( ( hist.y - fit.y )^2 / fit.y )
											})							
									ans			<- rbind(ans, data.table( height=ih, t.stat=t.stat, p= pchisq(t.stat, df=length(fit.x)-2, lower.tail=F), state=j, replicate=seq_along(t.stat), distr='poi', method='mo' ))						
									#	Geometric fit		
									try({
									mle			<- fitdistr(tmp, densfun="Geometric")
									fit.y		<- dgeom(fit.x, mle$estimate)
									if(plot) lines(fit.x, fit.y, col=cols[1], lty=ltys[3], lwd=lwds[3])
									t.stat		<- sum( ( hist.y - fit.y )^2 / fit.y )									
									ans			<- rbind(ans, data.table( height=ih, t.stat=t.stat, p= pchisq(t.stat, df=length(fit.x)-2, lower.tail=F), state=j, replicate=1, distr='geom', method='mle' ))
									}, silent=TRUE)
									t.stat		<- sapply(seq_along(pseudo.datasets), function(k)
											{
												fit.y		<- dgeom(fit.x, 1/(1+eMi[j,k]))
												if(plot) lines(fit.x, fit.y, col=util.fade.col(cols[4],0.4), lty=ltys[4], lwd=lwds[4])
												sum( ( hist.y - fit.y )^2 / fit.y )
											})									
									ans			<- rbind(ans, data.table( height=ih, t.stat=t.stat, p= pchisq(t.stat, df=length(fit.x)-2, lower.tail=F), state=j, replicate=seq_along(t.stat), distr='geom', method='mo' ))									
									#	NegBin fit		
									try({
									mle			<- fitdistr(tmp, densfun="negative binomial")
									fit.y		<- dnbinom(fit.x, size=mle$estimate['size'], mu=mle$estimate['mu'])
									if(plot) lines(fit.x, fit.y, col=cols[5], lty=ltys[5], lwd=lwds[5])
									t.stat		<- sum( ( hist.y - fit.y )^2 / fit.y )									
									ans			<- rbind(ans, data.table( height=ih, t.stat=t.stat, p= pchisq(t.stat, df=length(fit.x)-3, lower.tail=F), state=j, replicate=1, distr='negbin', method='mle' ))
									}, silent=TRUE)
									#	fit Neg Bin based on first two moments
									t.stat		<- sapply(seq_along(pseudo.datasets), function(k)
											{
												nbinom.p	<- 1 - eMi[j,k]/eMii[j,k]
												nbinom.r	<- eMi[j,k] * (1 - nbinom.p) / nbinom.p
												fit.y		<- dnbinom(fit.x, size=nbinom.r, prob=1-nbinom.p)
												if(plot) lines(fit.x, fit.y, col=util.fade.col(cols[6],0.4), lty=ltys[6], lwd=lwds[6])
												sum( ( hist.y - fit.y )^2 / fit.y )
											})									
									ans			<- rbind(ans, data.table( height=ih, t.stat=t.stat, p= pchisq(t.stat, df=length(fit.x)-3, lower.tail=F), state=j, replicate=seq_along(t.stat), distr='negbin', method='mo' ))
									if(plot && j==1)	legend('topright', bty='n', border=NA, lty=ltys, fill=cols, legend=legend.txt)
									ans
								})	
				)
				df.test	<- do.call('rbind',df.test)
				df.test				
			})	
	if(plot)	dev.off()		
	df.test	<- do.call('rbind',df.test)
	df.test			
}
################################################################################################
cg.sde.tiptypenumber.pseudolkl.for.param<- function(pseudo.datasets, heights, resume=1, file.r=NA, file.pdf=NA)
{
	require(data.table)
	if(resume && !is.na(file.r))
	{
		cat(paste("\ntry load ",file.r))
		options(show.error.messages = FALSE, warn=1)		
		readAttempt<-try(suppressWarnings(load(file.r)))						
		options(show.error.messages = TRUE)								
	}	
	if(!resume || inherits(readAttempt, "try-error"))
	{
		#collect cluster size distributions over replicates
		distr	<- lapply(seq_along(pseudo.datasets), function(k)
				{
					if(!k%%10) cat(paste('\nprocess pseudo.dataset #',k))
					pseudo.data	<- pseudo.datasets[[k]]				
					if(is.null(pseudo.data$distr))
					{
						bdt		<- pseudo.data$bdt					
						distr	<- calculate.cluster.size.distr.from.tree(bdt, heights)					
					}
					else
						distr	<- pseudo.data$distr
					cbind( distr, replicate=k )
				})
		distr	<- do.call("rbind",distr)	
		# plot distribution of cluster sizes for each tree height
		#cg.sde.tiptypenumber.pseudolkl.plotdistr(pseudo.datasets, distr, heights, file.pdf)	
		# explore potential pseudo likelihoods for each tree height
		df.test	<- cg.sde.tiptypenumber.pseudolkl.plotpotentiallkl(pseudo.datasets, distr, file.pdf)
		cat(paste("\nsave df.test to ",file.r))
		save(df.test, file=file.r)
	}
	df.test
}
################################################################################################
cg.sde.eM.pseudolkl.fit1DmMid.plot.pvalue<- function(df.fit, file)
{
	library(RColorBrewer)
	pseudolkls	<- unique( df.fit[, distr] )
	methods		<- unique( df.fit[, method] )
	parms.vary	<- unique( subset(df.fit, select=c(gi,bf)) )
	states		<- unique( df.fit[,state] )
	cols		<- brewer.pal(length(states), 'Set2')
	
	pdf(file=file, w=10, h=6)
	par(mfcol=c(length(pseudolkls),length(states)))
	par(mar=c(4,4,4,0.5))
	for(p in seq_len(nrow(parms.vary)))
	{
		cat(paste('\nprocess parameters #',p))
		i<- 2			
		dummy		<- sapply(seq_along(states), function(k)
				{
					dummy<- sapply(seq_along(pseudolkls), function(j)
							{							
								#plot moment p values							
								plot(1,1,type='n',bty='n',xlab='height',ylab='p-value',main=paste('(',parms.vary[p,gi],', ',parms.vary[p,bf],') ',pseudolkls[j],', ',methods[i], sep=''),ylim=c(0,1),xlim=range(df.fit[,height]))
								tmp			<- c(parms.vary[p,gi], parms.vary[p,bf]) 
								tmp			<- subset(df.fit, gi==tmp[1] & bf==tmp[2] &  distr==pseudolkls[j] & method==methods[i] & state==states[k])
								sapply(unique(tmp[,bdt]), function(r)
										{										
											tmp			<- subset(tmp, bdt==r)									
											points(tmp[,height],tmp[,p],col=util.fade.col(cols[k],0.3), pch=18)									
										})
								legend('topleft', legend=paste('state',states[k]),fill=cols[k],border=NA, bty='n')
								#add overall p values
								tmp			<- c(parms.vary[p,gi], parms.vary[p,bf])
								tmp			<- subset(df.fit, gi==tmp[1] & bf==tmp[2] &  distr==pseudolkls[j] & method==methods[1] & state==states[k])
								lines(tmp[,height],tmp[,p],pch=18)
							})	
				})						
	}
	dev.off()
}
################################################################################################
cg.sde.eM.pseudolkl.fit1DmMid<- function(df.obs, df.sim, file=NA)		
{
	require(nortest)
	heights			<- unique( df.obs[, height] )	
	state.cols.obs	<- which( grepl('state',colnames( df.obs )) ) 
	m				<- length( state.cols.obs )
	plot		<- !is.na(file)
	if(plot)
	{
		pdf(file=file, w=8, h=4)
		par(mfcol=c(1,m))
		cols		<- brewer.pal(3, 'Set1')[1:2]			
		ltys		<- c(1,1)
		lwds		<- rep(1.5,1.5)
		legend.txt	<- c('Normal mo','LogNormal mo')
	}		
	df.test<- lapply( seq_along(heights), function(ih)
			{
				df.test	<- lapply( seq_len(m), function(j)
						{
							su.obs	<- unclass( subset(df.obs, height==heights[ih])[,state.cols.obs[j],with=0] )[[1]]									
							su.sim	<- paste("subset(df.sim, height==heights[ih], select=c('bdt','replicate','state.",j,"') )",sep='')
							su.sim	<- eval( parse(text=su.sim) )
							setnames(su.sim,3,'state')
							hist.obs<- hist( su.obs, breaks=15, plot=FALSE )
							if(plot)	plot(hist.obs, freq=FALSE, main=paste('height=',round(heights[ih],d=1)), xlab=paste('eM',j), xlim=range(su.obs), ylim=c(0, max(hist.obs$density)*1.4), col='gray50', bty='n', border=NA)																										
							#	overall test for distribution of eMi							
							ans		<- data.table( height=heights[ih], state=j, bdt=1, distr='n', method='sf', p= sf.test(su.obs)$p.value )									
							ans		<- rbind(ans, data.table( height=heights[ih], state=j, bdt=1, distr='ln', method='sf', p= sf.test(log(su.obs))$p.value ))
							#	tests for distribution of eMi with moments specified by mMi 
							#	normal									
							tmp			<- su.sim[,{
										state.m		<- mean(state)
										state.sd	<- sd(state)												
										if(plot)	lines(hist.obs$mids, dnorm(hist.obs$mids, state.m, state.sd), col=util.fade.col(cols[1],0.25), lty=ltys[1], lwd=lwds[1])
										list(distr='n', method='mo-B',p= suppressWarnings(ks.test(su.obs, 'pnorm', mean=state.m, sd=state.sd)$p.value))	
									},by='bdt']
							ans		<- rbind(ans, cbind(height=heights[ih], state=j, tmp))				
							#	log normal									
							tmp			<- su.sim[,{
										state.m		<- mean(log(state))
										state.sd	<- sd(log(state))												
										if(plot)	lines(hist.obs$mids, dlnorm(hist.obs$mids, state.m, state.sd), col=util.fade.col(cols[2],0.25), lty=ltys[2], lwd=lwds[2])
										list(distr='ln', method='mo-B',p= suppressWarnings(ks.test(su.obs, 'plnorm', meanlog=state.m, sdlog=state.sd)$p.value))	
									},by='bdt']
							ans		<- rbind(ans, cbind(height=heights[ih], state=j, tmp))				
							ans
						})	
				if(plot)	legend('topright', bty='n', border=NA, lty=ltys, fill=cols, legend=legend.txt)
				do.call('rbind',df.test)					
			})
	df.test	<- do.call('rbind',df.test)
	if(plot)	dev.off()
	df.test
}
################################################################################################
cg.sde.eM.pseudolkl.plot.eMi.distr<- function(df.eMi, file=NA, lab='eM', method='2D')		
{
	if(!method%in%c('2D','1D'))	stop('unexpected method')
	heights		<- unique( df.eMi[, height] )	
	state.cols	<- which( grepl('state',colnames( df.eMi )) ) 
	m			<- length( state.cols )
	plot		<- !is.na(file)
	if(plot)
	{
		pdf(file=file, w=12, h=4)
		par(mfcol=c(1,m))
		cols		<- brewer.pal(3, 'Set1')
		cols		<- rep(cols,each=2)
		ltys		<- c(1,2,1,2,1,2)
		lwds		<- rep(1.5,6)
		legend.txt	<- c('Poisson mle','Poisson mo','Geometric mle','Geometric mo','NegBin mle', 'NegBin mo')
	}
	par(mfcol=c(1,m))
	dummy<- lapply( seq_along(heights), function(ih)
			{
				if(method=='1D')
				{
					dummy	<- sapply( seq_len(m), function(j)
							{
								tmp		<- unclass( subset(df.eMi, height==heights[ih])[,state.cols[j],with=0] )[[1]]
								hist.x	<- hist( tmp, main=paste('height=',round(heights[ih],d=1)), breaks=15, xlab=paste('eM',j), freq=FALSE, col='gray50', bty='n', border=NA, plot=plot )					
							})
				}
				if(method=='2D')
				{
					state.c	<- t( combn(3,2) ) 					
					dummy	<- sapply( seq_len(nrow(state.c)), function(j)
							{
								tmpx	<- unclass( subset(df.eMi, height==heights[ih])[,state.cols[state.c[j,1]],with=0] )[[1]]
								tmpy	<- unclass( subset(df.eMi, height==heights[ih])[,state.cols[state.c[j,2]],with=0] )[[1]]
								plot(tmpx, tmpy, main=paste('height=',round(heights[ih],d=1)), xlab=paste(lab,state.c[j,1]), ylab=paste(lab,state.c[j,2]), bty='n' )															
							})
				}
			})
	if(plot)	dev.off()
}
################################################################################################
cg.sde.eM.pseudolkl.plot.eMid.mMid.distr<- function(df.obs, df.sim, file,lab='delta eM')
{
	require(RColorBrewer)
	heights			<- unique( df.obs[, height] )	
	state.cols.obs	<- which( grepl('state',colnames( df.obs )) )
	state.cols.sim	<- which( grepl('state',colnames( df.sim )) )
	m				<- length( state.cols.obs )
	plot		<- !is.na(file)
	if(plot)
	{
		pdf(file=file, w=8, h=4)
		par(mfcol=c(1,m))
		cols		<- brewer.pal(3, 'Set1')[1]									
		legend.txt	<- c('calculated mo under theta0')
	}		
	state.c	<- t( combn(3,2) )
	dummy	<- lapply( seq_along(heights), function(ih)
			{					 					
				dummy	<- sapply( seq_len(nrow(state.c)), function(j)
						{
							tmpx	<- unclass( subset(df.sim, height==heights[ih])[,state.cols.sim[state.c[j,1]],with=0] )[[1]]
							tmpy	<- unclass( subset(df.sim, height==heights[ih])[,state.cols.sim[state.c[j,2]],with=0] )[[1]]
							plot(log(tmpx), log(tmpy), main=paste('height=',round(heights[ih],d=1)), xlab=paste('log',lab,state.c[j,1]), ylab=paste('log',lab,state.c[j,2]), bty='n', col=util.fade.col(cols[1],0.2), pch=18 )
							
							tmpx	<- unclass( subset(df.obs, height==heights[ih])[,state.cols.obs[state.c[j,1]],with=0] )[[1]]
							tmpy	<- unclass( subset(df.obs, height==heights[ih])[,state.cols.obs[state.c[j,2]],with=0] )[[1]]
							points(log(tmpx), log(tmpy) )																
						})
				legend('topleft', bty='n', border=NA, fill=cols, legend=legend.txt)								
			})
	if(plot) dev.off()
}
################################################################################################
cg.sde.eM.pseudolkl.for.param<- function(pseudo.datasets, resume=1, file.r=NA, file.pdf=NA)
{
	if(resume && !is.na(file.r))
	{
		cat(paste("\ntry load ",paste(file.r,'_fit1D.R',sep='')))
		options(show.error.messages = FALSE, warn=1)		
		readAttempt<-try(suppressWarnings(load(paste(file.r,'_fit1D.R',sep=''))))						
		options(show.error.messages = TRUE)								
	}	
	if(!resume || inherits(readAttempt, "try-error"))
	{
		library(RColorBrewer)
		
		m		<- nrow( pseudo.datasets[[1]]$eMi )
		heights	<- unique( pseudo.datasets[[1]]$mMi[,height] )
		
		#
		#	collate data structures
		#
		
		#collect empirical first moments over bdt replicates and for all heights
		df.eMi	<- lapply(seq_along(pseudo.datasets), function(k)
				{					
					tmp				<- pseudo.datasets[[k]]$eMi
					rownames(tmp)	<- paste('state',1:m,sep='.')					
					cbind( as.data.table(t( tmp )), height=heights, bdt=k )				
				})
		df.eMi	<- do.call('rbind', df.eMi)
		#collect empirical second moments over bdt replicates and for all heights
		df.eMij	<- lapply(seq_along(pseudo.datasets), function(k)
				{										
					tmp				<- sapply( seq_len( dim(pseudo.datasets[[k]]$eMij)[3] ), function(h)	pseudo.datasets[[k]]$eMij[,,h][upper.tri(pseudo.datasets[[k]]$eMij[,,1],diag=1)] )
					rownames(tmp)	<- sapply(1:m, function(i) paste('state.',i,1:m,sep=''))[upper.tri(pseudo.datasets[[k]]$eMij[,,1],diag=1)] 
					cbind( as.data.table(t( tmp )), height=heights, bdt=k )								
				})
		df.eMij	<- do.call('rbind', df.eMij)
		#collect calculated first moments under the true parameter theta0 over bdt replicates and 100 replicates for each bdt and for all heights
		df.mMi	<- lapply(seq_along(pseudo.datasets), function(k)
				{					
					cbind( pseudo.datasets[[k]]$mMi, bdt=k )										
				})
		df.mMi	<- do.call('rbind', df.mMi)
		#collect calculated second moments under the true parameter theta0 over bdt replicates and 100 replicates for each bdt and for all heights
		df.mMij	<- lapply(seq_along(pseudo.datasets), function(k)
				{					
					cbind( pseudo.datasets[[k]]$mMij, bdt=k )										
				})
		df.mMij	<- do.call('rbind', df.mMi)
		
		#
		#	clean up
		#
				
		#	clean up df.mMi: remove numerical inaccuracies in mMi
		cmd		<- paste('is.nan(state.',1:m,')',' | ','abs(state.',1:m,')<',1e-6,sep='', collapse=' | ')	
		cmd		<- paste('df.mMi[, which(',cmd,')]',sep='')
		tmp		<- eval( parse(text=cmd) )
		set(df.mMi, tmp, 'height', NA)				
		cmd		<- paste('{ ok.last= tail(which(!is.na(height)),1);    list(', paste('state.',1:m,'= state.',1:m,'[seq_len(ok.last)]',sep='', collapse=','),', height=height[seq_len(ok.last)]) }', sep='')
		cmd		<- paste('df.mMi[,',cmd,', by=c("bdt","replicate"),]')
		df.mMi	<- eval( parse(text=cmd) )
		
		#	remove negative changes in mMi -- ONLY FOR HOMOCHRONEOUS SAMPLING AND NONDECREASING POPSIZE
		cmd		<- paste('state.',1:m,'<=',0,sep='', collapse=' | ')	
		cmd		<- paste('df.mMi[, which(',cmd,')]',sep='')
		tmp		<- eval( parse(text=cmd) )
		set(df.mMi, tmp, 'height', NA)				
		cmd		<- paste('{ ok.last= tail(which(!is.na(height)),1);    list(', paste('state.',1:m,'= state.',1:m,'[seq_len(ok.last)]',sep='', collapse=','),', height=height[seq_len(ok.last)]) }', sep='')
		cmd		<- paste('df.mMi[,',cmd,', by=c("bdt","replicate"),]')
		df.mMi	<- eval( parse(text=cmd) )
		
		#
		#	get delta eMi and mMi
		#
		
		df.eMi[, dummy:= -height]
		setkey(df.eMi, bdt, dummy)
		cmd		<- paste('state.',1:m,'= -diff(state.',1:m,')',sep='', collapse=',')
		cmd		<- paste('df.eMi[, list(',cmd,', height=height[-1]), by="bdt"]',sep='')
		df.eMid	<- eval( parse(text=cmd) )
		setkey(df.eMid, bdt, height)
				
		df.mMi[, dummy:= -height]		
		setkey(df.mMi, bdt, replicate, dummy)
		cmd		<- paste('state.',1:m,'= -diff(state.',1:m,')',sep='', collapse=',')
		cmd		<- paste('df.mMi[, list(',cmd,', height=height[-1]), by=c("bdt","replicate")]',sep='')
		df.mMid	<- eval( parse(text=cmd) )
		setkey(df.mMid, bdt, replicate, height)

		if(0)	#	check if deltas make sense
		{
			ih<- 20
			subset( df.eMi, height==heights[ih])[, list(su1=summary(state.1), su2=summary(state.2), su3=summary(state.3), info=names(summary(state.1)))]		
			subset( df.mMi, height==heights[ih])[, list(su1=summary(state.1), su2=summary(state.2), su3=summary(state.3), info=names(summary(state.1)))]
			subset( df.eMid, height==heights[ih])[, list(su1=summary(state.1), su2=summary(state.2), su3=summary(state.3), info=names(summary(state.1)))]
			subset( df.mMid, height==heights[ih] & replicate==1)[, list(su1=summary(state.1), su2=summary(state.2), su3=summary(state.3), info=names(summary(state.1)))]
		}
		if(0)	#	a few diagnostic plots
		{
			cg.sde.eM.pseudolkl.plot.eMi.distr(df.eMi, lab='eM', method='2D', file=paste(file.pdf, '_eMi_2D.pdf',sep=''))
			cg.sde.eM.pseudolkl.plot.eMi.distr(df.eMid, lab='delta eM', method='1D', file=paste(file.pdf, '_eMid_1D.pdf',sep='') )
			cg.sde.eM.pseudolkl.plot.eMi.distr(df.eMid, lab='delta eM', method='2D', file=paste(file.pdf, '_eMid_2D.pdf',sep='') )
		}
		# compare different 1D pseudo lkls for each tree height
		df.fit1d	<- cg.sde.eM.pseudolkl.fit1DmMid(df.eMid, df.mMid, file=paste(file.pdf, '_eMid_1D_fit.pdf',sep=''))
					
		cg.sde.eM.pseudolkl.plot.eMid.mMid.distr(df.eMid, df.mMid, paste(file.pdf, '_eMid_2D.pdf',sep=''), lab='delta eM')	
		
		cat(paste("\nsave df.fit to ",paste(file.r,'_fit1D.R',sep='')))
		save(df.eMi, df.mMi, df.eMid, df.mMid, df.eMij, df.mMij, df.fit1d, file=paste(file.r,'_fit1D.R',sep=''))		
	}
	df.fit1d
}
################################################################################################
cg.sde.tiptypenumber.pseudolkl.for.param2<- function(pseudo.datasets, heights, resume=1, file.r=NA, file.pdf=NA)
{
	if(resume && !is.na(file.r))
	{
		cat(paste("\ntry load ",file.r))
		options(show.error.messages = FALSE, warn=1)		
		readAttempt<-try(suppressWarnings(load(file.r)))						
		options(show.error.messages = TRUE)								
	}	
	if(!resume || inherits(readAttempt, "try-error"))
	{
		#collect cluster size distributions over replicates
		distr	<- lapply(seq_along(pseudo.datasets), function(k)
				{
					if(!k%%10) cat(paste('\nprocess pseudo.dataset #',k))
					pseudo.data	<- pseudo.datasets[[k]]				
					if(is.null(pseudo.data$distr))
					{
						bdt		<- pseudo.data$bdt					
						distr	<- calculate.cluster.size.distr.from.tree(bdt, heights)					
					}
					else
						distr	<- pseudo.data$distr
					cbind( distr, replicate=k )
				})
		distr	<- do.call("rbind",distr)	
		# plot distribution of cluster sizes for each tree height
		#cg.sde.tiptypenumber.pseudolkl.plotdistr(pseudo.datasets, distr, heights, file.pdf)	
		# explore potential pseudo likelihoods for each tree height
		df.test	<- cg.sde.tiptypenumber.pseudolkl.plotpotentiallkl2(pseudo.datasets, distr, file=file.pdf)
		cat(paste("\nsave df.test to ",file.r))
		save(df.test, file=file.r)
	}
	df.test
}	
################################################################################################
cg.sde.tiptypenumber.pseudolkl.plot.pvalue<- function(df.test, file)
{
	library(RColorBrewer)
	pseudolkls	<- c('poi','geom','negbin')
	methods		<- c('mle','mo')
	states		<- unique(df.test[,state])
	cols		<- brewer.pal(length(states), 'Set2')
	
	
	pdf(file=file.pdf, w=10, h=6)
	par(mfcol=c(length(methods),length(pseudolkls)))
	par(mar=c(4,5,4,0.5))
	
	for(j in seq_along(pseudolkls))
		for(i in seq_along(methods))
		{
			plot(1,1,type='n',bty='n',xlab='height',ylab='chi2 pvalue of\ngoodness of fit in distribution',main=paste('distr=',pseudolkls[j],', method=',methods[i]),ylim=c(0,1),xlim=range(df.test[,height]))
			dummy		<- sapply(seq_along(states), function(k)
					{
						tmp			<- subset(df.test, distr==pseudolkls[j] & method==methods[i] & state==states[k])
						lines(tmp[,height],tmp[,p],col=cols[k])
						legend('bottomleft', legend=paste('state',states),fill=cols,border=NA, bty='n')
					})			
		}
	dev.off()
}
################################################################################################
cg.sde.tiptypenumber.pseudolkl.plot.pvalue2<- function(df.test, file)
{
	library(RColorBrewer)
	pseudolkls	<- c('poi','geom','negbin')
	methods		<- c('mle','mo')
	parms.vary	<- unique( subset(df.test, select=c(gi,bf)) )
	states		<- unique( df.test[,state] )
	cols		<- brewer.pal(length(states), 'Set2')
	
	
	pdf(file=file.pdf, w=10, h=6)
	par(mfcol=c(length(methods),length(pseudolkls)))
	par(mar=c(4,5,4,0.5))
	for(p in seq_len(nrow(parms.vary)))
	{
		cat(paste('\nprocess parameters #',p))
		for(j in seq_along(pseudolkls))
			for(i in seq_along(methods))
			{
				plot(1,1,type='n',bty='n',xlab='height',ylab='chi2 pvalue of\ngoodness of fit in distribution',main=paste(parms.vary[p,gi],',',parms.vary[p,bf],',',pseudolkls[j],',',methods[i]),ylim=c(0,1),xlim=range(df.test[,height]))
				dummy		<- sapply(seq_along(states), function(k)
						{
							tmp			<- c(parms.vary[p,gi], parms.vary[p,bf]) 
							tmp			<- subset(df.test, gi==tmp[1] & bf==tmp[2] &  distr==pseudolkls[j] & method==methods[i] & state==states[k])
							replicates	<- unique(tmp[,replicate])
							sapply(replicates, function(r)
									{
										tmp			<- subset(tmp, replicate==r)									
										lines(tmp[,height],tmp[,p],col=util.fade.col(cols[k],0.5))									
									})
						})			
				legend('bottomleft', legend=paste('state',states),fill=cols,border=NA, bty='n')
			}
	}
	dev.off()
}
################################################################################################
cg.sde.eM.pseudolkl.fit3D<- function()
{
	require(data.table)
	my.mkdir(HOME, 'MOMSDE' )
	dir.name		<- paste(HOME, 'MOMSDE',sep='/')
	S0				<- 5e3
	S0				<- 2e5
	phi				<- 0.5
	phi				<- 0.25
	#	parameter template
	parms.template 	<- list(	m=3, gamma0 = 1, gamma1 = 1/7, gamma2 = 1/2, mu = 1/30, b=.036, beta0 = 1+1/30, beta1=(1+1/30)/10, beta2=(1+1/30)/10, 
			S_0=S0, I0_0=1, I1_0=1, I2_0=1, alpha = 4, 
			times=seq(0, 50, by=.1), sampleTime=50, phi=phi)	
	INFECTEDNAMES 	<<- c('I0', 'I1', 'I2')
	run.old.version	<- 0
	resume			<- 1
	heights			<- seq(10, 50, length.out=50)
	parms.vary.i	<- 8
	#	bits of the model that are varied
	gi				<- c(0.5, 1, 3) 			#avg duration of I0
	bf				<- c(2,4,8,16) 				#dampening of beta0
	parms.vary		<- expand.grid(gi=gi, bf= bf)	
	parms.vary		<- cbind(parms.vary, gamma0= 1/parms.vary[,'gi'], beta1= parms.template[['beta0']]/parms.vary[,'bf'])
	#	get files to precomputed mM.sets	
	files.mMset		<- list.files(dir.name)
	files.mMset		<- files.mMset[ sapply( files.mMset, function(x)	grepl('pseudo.mM',x, fixed=1) & grepl(paste('S=',parms.template$S_0,sep=''),x, fixed=1) & grepl('set=',x, fixed=1) & grepl('.R',x, fixed=1)	) ]
	
	#	choose 'observed' data
	if(!is.na(parms.vary.i))
		parms.vary	<- parms.vary[parms.vary.i,,drop=0]
	
	df.fit		<- lapply(seq_len(nrow(parms.vary)),function(i)
			{				
				cat(paste('\nprocess',i,'gi=',parms.vary[i,gi],'bf=',parms.vary[i,bf]))
				# 	load pseudo data set
				file								<- paste(dir.name,'/',"pseudodata.S=",parms.template$S_0,"_gi=",parms.vary[i,'gi'],"_bf=",parms.vary[i,'bf'],".R",sep='')				
				if( is.na(file.info(file)$size) )	stop('cannot load data set')
				cat(paste("\nload file=",file))
				tmp									<- load(file)
				#	collect empirical first moments over bdt replicates and for all heights				
				m			<- nrow( pseudo.datasets[[1]]$eMi )
				heights		<- unique( pseudo.datasets[[1]]$mMi[,height] )								
				df.eMi		<- lapply(seq_along(pseudo.datasets), function(k)
								{					
									tmp				<- pseudo.datasets[[k]]$eMi
									rownames(tmp)	<- paste('state',1:m,sep='.')					
									cbind( as.data.table(t( tmp )), height=heights, bdt=k )				
								})
				df.eMi		<- do.call('rbind', df.eMi)
				#	get delta eMi and mMi				
				df.eMi[, dummy:= -height]
				setkey(df.eMi, bdt, dummy)
				cmd		<- paste('state.',1:m,'= -diff(state.',1:m,')',sep='', collapse=',')
				cmd		<- paste('df.eMi[, list(',cmd,', height=height[-1]), by="bdt"]',sep='')
				df.eMid	<- eval( parse(text=cmd) )
				setkey(df.eMid, bdt, height)
				
				
				
				
				
				rownames(tmp)	<- sapply(1:m, function(i) paste('state.',i,1:m,sep=''))[upper.tri(diag(m),diag=1)]
				
				subset( tmp, height== unique(tmp[,height])[26])	
				
				cov( tmp_h )
				
			})
	
}
################################################################################################
cg.sde.eM.pseudolkl<- function()
{
	my.mkdir(HOME, 'MOMSDE' )
	dir.name		<- paste(HOME, 'MOMSDE',sep='/')
	S0				<- 5e3
	S0				<- 2e5
	phi				<- 0.5
	phi				<- 0.25
	#parameter template
	parms.template 	<- list(	m=3, gamma0 = 1, gamma1 = 1/7, gamma2 = 1/2, mu = 1/30, b=.036, beta0 = 1+1/30, beta1=(1+1/30)/10, beta2=(1+1/30)/10, 
			S_0=S0, I0_0=1, I1_0=1, I2_0=1, alpha = 4, 
			times=seq(0, 50, by=.1), sampleTime=50, phi=phi)	
	INFECTEDNAMES 	<<- c('I0', 'I1', 'I2')
	run.old.version	<- 0
	resume			<- 0
	heights			<- seq(10, 50, length.out=50)
	
	if(resume)
	{
		file.r	<- 	paste(dir.name,'/',"pseudodistr.potential.lkl.S=",parms.template$S_0,"_fit1D.R",sep='')
		cat(paste("\ntry load ",file.r))
		options(show.error.messages = FALSE, warn=1)		
		readAttempt<-try(suppressWarnings(load(file.r)))						
		options(show.error.messages = TRUE)								
	}	
	if(!resume || inherits(readAttempt, "try-error"))
	{
		#bits of the model that are varied
		gi			<- c(0.5, 1, 3) 			#avg duration of I0
		bf			<- c(2,4,8,16) 				#dampening of beta0
		parms.vary	<- expand.grid(gi=gi, bf= bf)	
		parms.vary	<- cbind(parms.vary, gamma0= 1/parms.vary[,'gi'], beta1= parms.template[['beta0']]/parms.vary[,'bf'])
		
		df.fit		<- lapply(seq_len(nrow(parms.vary)),function(i)
				{
				
					cat(paste('\nprocess',i))
					# load pseudo data set
					parms								<- parms.template
					parms[c('gamma0','beta1','beta2')]	<- parms.vary[i,c('gamma0','beta1','beta1')]	
					file								<- paste("pseudodata.S=",parms.template$S_0,"_gi=",parms.vary[i,'gi'],"_bf=",parms.vary[i,'bf'],".R",sep='')
					file								<- paste(dir.name,file,sep='/')
					if( is.na(file.info(file)$size) )
						return(data.table(height=numeric(0), t.stat=numeric(0), p=numeric(0), state=numeric(0), replicate=numeric(0), distr=character(0), method=character(0), gi=numeric(0), bf=numeric(0)))
					cat(paste("\nload file=",file))
					tmp									<- load(file)
					# get pseudo.lkl p-values of potential pseudo lkl densities to the empirical density
					
					file.r			<- paste(dir.name,'/',"pseudo.eM.lkl.S=",parms.template$S_0,"_gi=",parms.vary[i,'gi'],"_bf=",parms.vary[i,'bf'],sep='')								
					file.pdf		<- paste(dir.name,'/',"pseudo.eM.lkl.S=",parms.template$S_0,"_gi=",parms.vary[i,'gi'],"_bf=",parms.vary[i,'bf'],sep='')				
					df.fit			<- cg.sde.eM.pseudolkl.for.param(pseudo.datasets, resume=0, file.r=file.r, file.pdf=file.pdf)
					df.fit[, gi:=parms.vary[i,'gi']]
					df.fit[, bf:=parms.vary[i,'bf']]
					df.fit
				})
		df.fit	<- do.call('rbind', df.fit)
		file.r	<- 	paste(dir.name,'/',"pseudo.eM.lkl.S=",parms.template$S_0,"_fit1D.R",sep='')
		cat(paste('\nsave df.fit to',file.r))
		save(df.fit, file=file.r)	
	}
		
	
	file.pdf	<- 	paste(dir.name,'/',"pseudo.eM.lkl.S=",parms.template$S_0,"_fit1Dpvalue.pdf",sep='')
	cg.sde.eM.pseudolkl.fit1DmMid.plot.pvalue(df.fit, file.pdf)	

}
################################################################################################
cg.sde.eval.pseudolkl.to.empirical.clustersize.distribution<- function()
{
	my.mkdir(HOME, 'MOMSDE' )
	dir.name		<- paste(HOME, 'MOMSDE',sep='/')
	S0				<- 5e3
	S0				<- 2e5
	phi				<- 0.5
	phi				<- 0.25
	#parameter template
	parms.template 	<- list(	m=3, gamma0 = 1, gamma1 = 1/7, gamma2 = 1/2, mu = 1/30, b=.036, beta0 = 1+1/30, beta1=(1+1/30)/10, beta2=(1+1/30)/10, 
								S_0=S0, I0_0=1, I1_0=1, I2_0=1, alpha = 4, 
								times=seq(0, 50, by=.1), sampleTime=50, phi=phi)	
	INFECTEDNAMES 	<<- c('I0', 'I1', 'I2')
	run.old.version	<- 0
	resume			<- 0
	heights			<- seq(10, 50, length.out=50)
	
	if(resume)
	{
		file.r	<- 	paste(dir.name,'/',"pseudodistr.potential.lkl.S=",parms.template$S_0,"_dftest2.R",sep='')
		cat(paste("\ntry load ",file.r))
		options(show.error.messages = FALSE, warn=1)		
		readAttempt<-try(suppressWarnings(load(file.r)))						
		options(show.error.messages = TRUE)								
	}	
	if(!resume || inherits(readAttempt, "try-error"))
	{
		#bits of the model that are varied
		gi			<- c(0.5, 1, 3) 			#avg duration of I0
		bf			<- c(2,4,8,16) 				#dampening of beta0
		parms.vary	<- expand.grid(gi=gi, bf= bf)	
		parms.vary	<- cbind(parms.vary, gamma0= 1/parms.vary[,'gi'], beta1= parms.template[['beta0']]/parms.vary[,'bf'])
		
		df.test		<- lapply(seq_len(nrow(parms.vary)),function(i)
				{
					cat(paste('\nprocess',i))
					# load pseudo data set
					parms								<- parms.template
					parms[c('gamma0','beta1','beta2')]	<- parms.vary[i,c('gamma0','beta1','beta1')]	
					file								<- paste("pseudodata.S=",parms.template$S_0,"_gi=",parms.vary[i,'gi'],"_bf=",parms.vary[i,'bf'],".R",sep='')
					file								<- paste(dir.name,file,sep='/')
					if( is.na(file.info(file)$size) )
						return(data.table(height=numeric(0), t.stat=numeric(0), p=numeric(0), state=numeric(0), replicate=numeric(0), distr=character(0), method=character(0), gi=numeric(0), bf=numeric(0)))
					cat(paste("\nload file=",file))
					tmp									<- load(file)
					# get pseudo.lkl p-values of potential pseudo lkl densities to the empirical density
					if(run.old.version)	
					{
						#old code fitting the pseudo lkl based on the pooled eMi.avg and eMii.avg
						file.r		<- paste(dir.name,'/',"pseudodistr.potential.lkl.S=",parms.template$S_0,"_gi=",parms.vary[i,'gi'],"_bf=",parms.vary[i,'bf'],"_dftest.R",sep='')
						file.pdf	<- paste(dir.name,'/',"pseudodistr.potential.lkl.S=",parms.template$S_0,"_gi=",parms.vary[i,'gi'],"_bf=",parms.vary[i,'bf'],"_distr.pdf",sep='')
						df.test		<- cg.sde.tiptypenumber.pseudolkl.for.param(pseudo.datasets, heights, resume=1, file.r=file.r, file.pdf=file.pdf)
					}
					else
					{
						#fit the pseudo lkl based on a single eMi and eMii
						file.r		<- paste(dir.name,'/',"pseudodistr.potential.lkl.S=",parms.template$S_0,"_gi=",parms.vary[i,'gi'],"_bf=",parms.vary[i,'bf'],"_dftest2.R",sep='')
						file.pdf	<- paste(dir.name,'/',"pseudodistr.potential.lkl.S=",parms.template$S_0,"_gi=",parms.vary[i,'gi'],"_bf=",parms.vary[i,'bf'],"_distr2.pdf",sep='')
						df.test		<- cg.sde.tiptypenumber.pseudolkl.for.param2(pseudo.datasets, heights, resume=0, file.r=file.r, file.pdf=file.pdf)
					}
					df.test[, gi:=parms.vary[i,'gi']]
					df.test[, bf:=parms.vary[i,'bf']]
					df.test
				})
		df.test	<- do.call('rbind', df.test)
		file.r	<- 	paste(dir.name,'/',"pseudodistr.potential.lkl.S=",parms.template$S_0,"_dftest2.R",sep='')
		cat(paste('\nsave df.tets to',file.r))
		save(df.test, file=file.r)
	}
	if(run.old.version)
	{
		file.pdf	<- paste(dir.name,'/',"pseudodistr.potential.lkl.S=",parms.template$S_0,"_pvalues.pdf",sep='')
		cg.sde.tiptypenumber.pseudolkl.plot.pvalue(df.test, file=file.pdf)
	}
	else
	{
		file.pdf	<- paste(dir.name,'/',"pseudodistr.potential.lkl.S=",parms.template$S_0,"_pvalues2.pdf",sep='')
		cg.sde.tiptypenumber.pseudolkl.plot.pvalue2(df.test, file=file.pdf)
		
	}
	
	
	
	
	
	#	for every single bdt, the cluster size distribution breaks down deep in the tree, need h<40
	#	for every single bdt, the variance is hard to estimate close to the tips, need h>35 or more tips
}
################################################################################################
cg.sde.get.mM.collect.runs<- function()
{
	require(data.table)
	parms.template 	<- list(	m=3, gamma0 = 1, gamma1 = 1/7, gamma2 = 1/2, mu = 1/30, b=.036, beta0 = 1+1/30, beta1=(1+1/30)/10, beta2=(1+1/30)/10, 
								S_0=2e5, I0_0=1, I1_0=1, I2_0=1, alpha = 4, 
								times=seq(0, 50, by=.1), sampleTime=50, phi=0.25, mM.replicate= 100)	
		
	dir.name		<- paste(HOME, '/MOMSDE', '/pseudo.mM.S=', parms.template$S_0,sep='')	
	files			<- list.files(dir.name)	
	files			<- files[ sapply( files, function(x)	grepl('gi=',x, fixed=1)	) ]	
	mM.nsets		<- ceiling( length(files)/1e3 )
	mM.files		<- data.table(file=files, mM.set=ceiling(seq_along(files)/1e3))
	m				<- parms.template$m
	
	#	collect simulated mM for various parms 
	dummy	<- sapply(unique(mM.files[,mM.set]), function(mM.set.i)
			{
				cat(paste('\nprocess set',mM.set.i))
				mM.files.set	<- subset(mM.files, mM.set==mM.set.i)
				mM.set			<- lapply( seq_len(nrow(mM.files.set)), function(i)
						{							
							load( paste(dir.name, mM.files.set[i,file], sep='/') )
							mMid	<- sim$mMid
							#	fix cleaning of mMid
							#	clean up df.mMi: remove numerical inaccuracies in mMi
							cmd		<- paste('state.',1:m,'<',0,sep='', collapse=' | ')	
							cmd		<- paste('mMid[, which(',cmd,')]',sep='')
							tmp		<- eval( parse(text=cmd) )
							set(mMid, tmp, 'height', NA)	
							if(mMid[, any( is.na(height))])
							{
								cmd		<- paste('{ ok.last<- which(is.na(height));  ok.last<- ifelse(length(ok.last), ok.last-1, length(height));  list(', paste('state.',1:m,'= state.',1:m,'[seq_len(ok.last)]',sep='', collapse=','),', height=height[seq_len(ok.last)]) }', sep='')
								cmd		<- paste('mMid[,',cmd,', by="replicate"]')
								mMid	<- eval( parse(text=cmd) )
							}
							#
							mMid[, gi:= 1/sim$parms$gamma0]
							mMid[, bf:= sim$parms$beta0 / sim$parms$beta1]
							mMid				
						})
				mM.set			<- do.call('rbind', mM.set)
				file			<- paste(HOME, '/MOMSDE', '/pseudo.mM.S=', parms.template$S_0,'.set=',mM.set.i,'.R',sep='')
				save(mM.set, file=file)				
			})
	
	#	clean up: keep only parms for which we have >50 replicates for heights in 10-40
	dir.name		<- paste(HOME, 'MOMSDE', sep='/')
	files			<- list.files(dir.name)
	files			<- files[ sapply( files, function(x)	grepl('pseudo.mM',x, fixed=1) & grepl(paste('S=',parms.template$S_0,sep=''),x, fixed=1) & grepl('set=',x, fixed=1) & grepl('.R',x, fixed=1)	) ]	
	dummy			<- lapply( seq_along(files), function(mM.set.i)
			{
				file			<- paste(dir.name,files[mM.set.i],sep='/')
				load(file)
				
				parms.vary		<- unique( subset( mM.set, select=c(gi, bf) ) )
				
				mM.set.select	<- mM.set[, list( height.ok.n=length(which( height<40 & height>10 ) )), by=c('gi', 'bf', 'replicate')]
				mM.set.select	<- mM.set.select[, list( replicate.n.with.height.ok=length(which(height.ok.n>=25)) ), by=c('gi', 'bf')]
				mM.set.select	<- subset( mM.set.select, replicate.n.with.height.ok>50 )
				#
				cat(paste('\nplot mM.set.select to file=',paste(substr(file, 1, nchar(file)-2),'_mMselect.pdf', sep='')))
				pdf( paste(substr(file, 1, nchar(file)-2),'_mMselect.pdf', sep=''), 4, 4)
				plot( mM.set.select[,gi], mM.set.select[, bf], bty='n', pch=18, xlab='gi', ylab='bf')
				dev.off()
				#
				mM.set			<- merge( subset( mM.set.select, select=c(gi, bf) ), mM.set, by=c('gi', 'bf'))
				cat(paste('\nsave cleaned mM.set to file', file))
				save(mM.set, file=file)
			})
	
	#	construct log normal pseudo likelihoods
	dir.name		<- paste(HOME, 'MOMSDE', sep='/')
	files			<- list.files(dir.name)
	files			<- files[ sapply( files, function(x)	grepl('pseudo.mM',x, fixed=1) & grepl(paste('S=',parms.template$S_0,sep=''),x, fixed=1) & grepl('set=',x, fixed=1) & grepl('.R',x, fixed=1)	) ]	
	dummy			<- lapply( seq_along(files), function(mMset.i)
			{
				#mMset.i	<- 1
				file	<- paste(dir.name, files[mMset.i], sep='/')
				cat(paste('\nprocess',file))
				load( file  )
				#	log mMid	
				mM.states.col	<- which( grepl('state',colnames( mM.set )) )
				for(j in mM.states.col)
					set(mM.set, NULL, colnames(mM.set)[j], log( unclass( mM.set[, j, with=0] )[[1]]) )				
				#	for each parms and each height, construct lognormal pseudolikelihood parameters
				setkey( mM.set, gi, bf, height)
				parms	<- unique( subset( mM.set, select=c(gi, bf)))
				mM.ln	<- lapply(seq_len(nrow(parms)), function(parms.j)
						{
							tmp_p	<- subset(mM.set, gi==parms[parms.j,gi] & bf==parms[parms.j,bf])				
							cmd		<- paste('{ 	tmp_h		<- as.matrix(cbind(',paste('state.',1:m, collapse=', ', sep=''),')); 
											ans			<- as.vector(cov( tmp_h )[ upper.tri(diag(nrow=m, ncol=m), diag=TRUE) ]); 
											names(ans)	<- sapply(1:m, function(i) paste("cov.",i,1:m,sep=""))[upper.tri(diag(m),diag=1)];
											as.list( c( colMeans( tmp_h ), ans) ) }',sep='')
							cmd		<- paste('tmp_p[,',cmd,', by="height"]', sep='')
							mM.ln_p	<- eval( parse(text=cmd) )
							mM.ln_p[, gi:= parms[parms.j,gi] ]
							mM.ln_p[, bf:= parms[parms.j,bf] ]
							mM.ln_p
						})
				mM.ln	<- do.call('rbind', mM.ln)	
				file	<- paste(dir.name,'/',sub('mM','LNmM',files[mMset.i]),sep='')
				cat(paste('\nsave mM.ln to file=',file))
				save(mM.ln, file=file)
			})
	
	#check stuff
	tmp		<- subset(mM.set, gi==0.5005 & bf==2.42)
	tmp_h	<- subset( tmp, height==unique( tmp[,height] )[1] )
}
################################################################################################
cg.sde.get.mM<- function()
{	
	index		<- NA
	if(exists("argv"))
	{
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,2),
									i= return(as.numeric(substr(arg,4,nchar(arg)))),NA)	}))
		if(length(tmp)>0) index<- tmp[1]						
	}	
	
	my.mkdir(HOME, 'MOMSDE' )
	dir.name		<- paste(HOME, 'MOMSDE',sep='/')	
	n.parms			<- 2e4
	resume			<- 1
	parms.template 	<- list(	m=3, gamma0 = 1, gamma1 = 1/7, gamma2 = 1/2, mu = 1/30, b=.036, beta0 = 1+1/30, beta1=(1+1/30)/10, beta2=(1+1/30)/10, 
								S_0=2e5, I0_0=1, I1_0=1, I2_0=1, alpha = 4, 
								times=seq(0, 50, by=.1), sampleTime=50, phi=0.25, mM.replicate= 100)	
	INFECTEDNAMES 	<<- c('I0', 'I1', 'I2')
	file.parms		<- paste(dir.name, '/pseudo.mM.S=',parms.template$S_0,'_parms.R', sep='')
	file.mM			<- paste(dir.name, '/pseudo.mM.S=',parms.template$S_0,'_mM.R', sep='')
	
	#
	#	create parms.vary if necessary
	#
	if(resume)
	{		
		cat(paste("\ntry load ",file.parms))
		options(show.error.messages = FALSE, warn=1)		
		readAttempt	<- try(suppressWarnings(load(file.parms)))						
		options(show.error.messages = TRUE)								
	}	
	if(!resume || inherits(readAttempt, "try-error"))	
	{
		gi			<- round(runif(n.parms, 0.2, 4), d=4) 				#avg duration of I0
		bf			<- round(runif(n.parms, 1,16), d=2) 				#dampening of beta0
		parms.vary	<- unique( cbind(gi=gi, bf= bf), MARGIN=1 )	
		parms.vary	<- cbind(parms.vary, gamma0= 1/parms.vary[,'gi'], beta1= parms.template[['beta0']]/parms.vary[,'bf'])
		save(parms.vary, file=file.parms)
	}
		
	#	if moments are calculated in parallel, reduce parms.vary	
	if(!is.na(index))
	{
		parms.vary	<- parms.vary[index,,drop=0]
	}
	#	calculate moments
	sims	<- lapply(seq_len(nrow(parms.vary)), function(i)
			{
				cat(paste("\nprocess gi=",parms.vary[i,'gi'],", bf=",parms.vary[i,'bf'],', i=',i))
				parms.template[c('gamma0','beta1','beta2')]	<- parms.vary[i,c('gamma0','beta1','beta1')]	
				parms			<- parms.template
				file			<- paste(dir.name,'/',"pseudo.mM.S=",parms.template$S_0,"_gi=",parms.vary[i,'gi'],"_bf=",parms.vary[i,'bf'],".R",sep='')				
				sim				<- cg.sde.get.mM.for.param(parms, file=file)				
			})	
	if(is.na(index))
	{
		cat(paste("\nsave pseudo data sets to file=",file.mM))
		save(sims, file=file.mM)		
	}
}
################################################################################################
cg.sde.nsim.mM<- function()
{	
	my.mkdir(HOME, 'MOMSDE' )
	dir.name		<- paste(HOME, 'MOMSDE',sep='/')
	
	#~ parms_truth <<- list(gamma0 = 1, gamma1 = 1/7, gamma2 = 1/2, mu = 1/30, b=.036, beta0 = 0.775, beta1=0.08, beta2=0.08, S0=2500, alpha = .05) 
	parms_truth 	<<- list(	m=3, gamma0 = 1, gamma1 = 1/7, gamma2 = 1/2, mu = 1/30, b=.036, beta0 = 1+1/30, beta1=(1+1/30)/10, beta2=(1+1/30)/10, 
								S_0=2e5, I0_0=1, I1_0=1, I2_0=1, alpha = 4,
								times=seq(0, 50, by=.1), sampleTime=50, phi=0.25) 
	FGYPARMS 		<<- parms_truth
	INFECTEDNAMES 	<<- c('I0', 'I1', 'I2')
	nsims			<- 1
	
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
	file			<- paste('nsim.mM.epidemic.S=',parms_truth$S_0,'.pdf',sep='')
	file			<- paste(dir.name,file,sep='/')
	dummy			<- solve.model.set.fgy(parms_truth, file) 
	Y.sampleTime 	<- Y.(parms_truth$sampleTime)
	m				<- 3
	stateIndices 	<- rep( 1:m, round( parms_truth$phi * Y.sampleTime ) ) # sample each of three states in proportion to size
	n 				<- length(stateIndices)
	sampleTimes 	<- rep(parms_truth$sampleTime, n)
	sampleStates 	<- diag(m)[stateIndices,]
	cat(paste("\nsimulate structured coalescent tree, including tip states (0/1)"))
	coalescentTree_time <- system.time({
				bdt <- simulatedBinaryDatedTree(parms_truth$sampleTime, sampleStates, discretizeRates=TRUE) 
			})
	print(coalescentTree_time)
	# calculate empirical stats
	heights 		<- seq(0, 50, length.out=50)
	cat(paste("\ncalculate moments of observed tree by averaging over lineages; tip states 0/1"))
	eM 				<- calculate.cluster.size.moments.from.tree(bdt, heights)
	eM				<- eM$Mij
	cat(paste("\ncalculate moments under model"))
	mM				<- cg.sde.modelMoments(parms_truth, parms_truth$sampleTime, sampleStates, solve.model.set.fgy, nsims)
	mM				<- mM$mMij
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
