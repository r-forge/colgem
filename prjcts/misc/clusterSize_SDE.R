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
	mM <- list()
	#for (isim in 1:nsims)
	#{
	#	mM.time <- system.time({  
	#				modelMoments = mM_i <- calculate.cluster.size.moments.from.model(sampleTime, sampleStates , timeResolution = 50, discretizeRates=TRUE, fgyResolution = 100 , integrationMethod = 'rk4')
	#			})
	#}
	mM <- lapply(1:nsims, function(isim) 
			{
				solve.model.set.fgy(parms)
				#timeResolution = 50; discretizeRates=TRUE; fgyResolution = 100; integrationMethod = 'rk4'
				calculate.cluster.size.moments.from.model(sampleTime, sampleStates , timeResolution = 50, discretizeRates=TRUE, fgyResolution = 100 , integrationMethod = 'rk4');
			})
	list(mM=mM)	
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
	cg.sde.eval.pseudolkl.to.empirical.clustersize.distribution()
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
cg.sde.get.pseudodata.for.param<- function(parms, bdt.n= 1e2, bdt.heights=seq(0, 50, length.out=50), distr.heights=seq(10, 50, length.out=50))
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
	
	#simulated btd.n pseudo data sets
	pseudo.datasets	<- lapply(seq_len(bdt.n), function(i)
			{
				pseudo.data					<- solve.model.set.fgy(parms_truth)	
				pseudo.data$stateIndices 	<- rep( 1:m, round( parms$phi * pseudo.data$Y( parms$sampleTime ) ) ) # sample each of three states in proportion to size	 
				pseudo.data$sampleTimes 	<- rep(parms$sampleTime, length(pseudo.data$stateIndices) )
				pseudo.data$sampleStates 	<- diag(m)[pseudo.data$stateIndices,]
				pseudo.data$bdt 			<- simulatedBinaryDatedTree(parms$sampleTime, pseudo.data$sampleStates, discretizeRates=TRUE)				
				pseudo.data$distr			<- calculate.cluster.size.distr.from.tree(pseudo.data$bdt, distr.heights)
				tmp			 				<- calculate.cluster.size.moments.from.tree(pseudo.data$bdt, bdt.heights)				
				pseudo.data$eMi 			<- tmp$Mi
				pseudo.data$eMij 			<- tmp$Mij
				pseudo.data
			})
	pseudo.datasets	
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
cg.sde.get.pseudodata<- function()
{	
	index	<- NA
	if(exists("argv"))
	{
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,2),
									i= return(as.numeric(substr(arg,4,nchar(arg)))),NA)	}))
		if(length(tmp)>0) index<- tmp[1]
	}	
	
	my.mkdir(HOME, 'MOMSDE' )
	dir.name		<- paste(HOME, 'MOMSDE',sep='/')	
	#parameter template
	parms.template 	<- list(	m=3, gamma0 = 1, gamma1 = 1/7, gamma2 = 1/2, mu = 1/30, b=.036, beta0 = 1+1/30, beta1=(1+1/30)/10, beta2=(1+1/30)/10, 
								S_0=2e5, I0_0=1, I1_0=1, I2_0=1, alpha = 4, 
								times=seq(0, 50, by=.1), sampleTime=50, phi=0.25)	
	INFECTEDNAMES 	<<- c('I0', 'I1', 'I2')
	
	#bits of the model that are varied
	gi			<- c(0.5, 1, 3) 			#avg duration of I0
	bf			<- c(2,4,8,16) 				#dampening of beta0
	parms.vary	<- expand.grid(gi=gi, bf= bf)	
	parms.vary	<- cbind(parms.vary, gamma0= 1/parms.vary[,'gi'], beta1= parms.template[['beta0']]/parms.vary[,'bf'])
	
	if(!is.na(index))
	{
		parms.vary	<- parms.vary[index,,drop=0]
	}
	#generate pseudo data sets
	dummy	<- lapply(seq_len(nrow(parms.vary)), function(i)
			{
				cat(paste("\nprocess gi=",parms.vary[i,'gi'],"bf=",parms.vary[i,'bf']))
				parms.template[c('gamma0','beta1','beta2')]	<- parms.vary[i,c('gamma0','beta1','beta1')]	
				parms			<-  parms.template
				pseudo.datasets	<- cg.sde.get.pseudodata.for.param(parms, bdt.n= 2, bdt.heights=seq(0, 50, length.out=50))
				file			<- paste(dir.name,'/',"pseudodata.S=",parms.template$S_0,"_gi=",parms.vary[i,'gi'],"_bf=",parms.vary[i,'bf'],".R",sep='')				
				cat(paste("\nsave pseudo data sets to file=",file))
				save(pseudo.datasets, parms, file=file)				
			})			
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
