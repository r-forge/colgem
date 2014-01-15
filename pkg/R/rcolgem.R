#' this file contains all R functions of the rcolgem package
#' expects F.(t), G.(t), and Y.(t) to be in namespace
#' @import ape
#' @import deSolve
#' @import pomp

#~ TODO finite size correction for pair (i,j) of lineages at internal node (see written notes): 
#~	if i transmits and is type k: 
#~ 		pjk -> pjk * (Yk-1)/Yk
#~ 		pjl (l\neq k) -> pjl * Yk/(Yk-pjk)
#~ TODO option to correct for direct ancestor sampling if doing serial samples and there is a 0-branch length
#~ 	add these terms to the likelihood:
#~ 		if 0 bl at s_i: \sum_k pik Ak / Yk
#~ 		else: \sum_k pik (1-Ak/Yk)
#~ TODO validate input & raise warnings
#~ 	check sampleTimes compatible with edge.length; FGY functions defined over length of tree;
#~ TODO adapt coalescent simulator to heterochronous sample
#~ TODO add option to switch to semi-structured coalescent if difference in line states below threshold
#~ 			TODO how to quantify difference in line state vectors
#~ 				max or mean euclidean distance to centroid of points on 1-simplex? 
#~ 			solve deterministic Ak to max time, likelihood based on total rather than pairwise coalescent rate
#~ TODO test heterochronous samples
#~ TODO a rare bug : Error in if (sum(tree$mstates[i, ]) <= 0) { : missing value where TRUE/FALSE needed

#require(ape)
#require(deSolve)

#' @export
SIMULATIONTIMERESOLUTION<- 1e+04

#PRELIMINARIES 
.calculate.heights <- function(phylo){
	phylo$maxSampleTime  = phylo$maxSampleTime <- max(phylo$sampleTimes)
	heights <- rep(0, (phylo$Nnode + length(phylo$tip.label)) )
	heights[1:length(phylo$sampleTimes)] <- phylo$maxSampleTime - phylo$sampleTimes
	curgen <- 1:length(phylo$sampleTimes)
	while( length(curgen) > 0) { 
		nextgen <- c()
		icurgenedges <- which(  phylo$edge[,2] %in% curgen  )
		for (i in icurgenedges){
			u<- phylo$edge[i,1]
			v<- phylo$edge[i,2]
			heights[u] <- phylo$edge.length[i] + heights[v]
			#if (is.na(heights[u])) browser()
			nextgen <- c(nextgen, u)
		}
		curgen <- unique(nextgen)
	}
	phylo$heights <- heights
	return(phylo)
}

.calculate.edgemap <- function(phylo){
	inEdgeMap <- rep(-1, length((phylo$Nnode + length(phylo$tip.label))))
	outEdgeMap <- matrix(-1, (phylo$Nnode + length(phylo$tip.label)), 2)
	parent <- 1:(phylo$Nnode + length(phylo$tip.label)) 
	daughters <- matrix(-1, (phylo$Nnode + length(phylo$tip.label)), 2)
	for (u in 1:(phylo$Nnode + length(phylo$tip.label))){
		#if (u!=length(phylo$tip.label)+1){ #if u not root
		tryCatch({ inEdgeMap[u] <- which( phylo$edge[,2]==u ) }, error = function(e) {inEdgeMap[u] <- u} )
		#} else{ 
		#	inEdgeMap[u] <- u
		#}
		if (u > length(phylo$tip.label)){
			outEdgeMap[u,] <- which( phylo$edge[,1]==u ) 
			daughters[u,] <- phylo$edge[outEdgeMap[u,],2]
		} else{ 
			outEdgeMap[u,] <- c(NA, NA)
			daughters[u,] <- c(NA, NA)
		}
		parent[u] <- phylo$edge[inEdgeMap[u],1]
	}
	phylo$inEdgeMap <- inEdgeMap
	phylo$outEdgeMap <- outEdgeMap
	phylo$parent = phylo$parents <- parent
	phylo$daughter = phylo$daughters <- daughters
	phylo$parentheight = phylo$parentheights <- phylo$heights[parent[1:(phylo$Nnode + length(phylo$tip.label))]]
	return(phylo)
}

.initialize.states <- function(phylo)
{
	phylo$m =m <- dim(phylo$sampleStates)[2]
	phylo$lstates <- matrix(-1, (phylo$Nnode + length(phylo$tip.label)), m)
	phylo$mstates <- matrix(-1, (phylo$Nnode + length(phylo$tip.label)), m)
	phylo$ustates <- matrix(-1, (phylo$Nnode + length(phylo$tip.label)), m)
	phylo$lstates[1:nrow(phylo$sampleStates),] <- phylo$sampleStates
	phylo$mstates[1:nrow(phylo$sampleStates),] <- phylo$sampleStates
	phylo$ustates[1:nrow(phylo$sampleStates),] <- phylo$sampleStates
	return(phylo)
}


#' Create binary dated tree
#' binaryDatedTree class, includes heights for each node and other helper variables
#' @export
binaryDatedTree <- function( x, ...) UseMethod("binaryDatedTree")

binaryDatedTree.default <- function( phylo, sampleTimes, sampleStates){
	if (phylo$Nnode != length(phylo$tip.label) - 1 ) { stop('Object class phylo is not a binary tree.') }
	phylo$sampleTimes <- sampleTimes
	phylo$sampleStates <- sampleStates
	phylo <- .calculate.heights(phylo)
	phylo <- .calculate.edgemap(phylo)
	phylo <- .initialize.states(phylo)
	phylo$coalescentRates <- rep(0, (phylo$Nnode + length(phylo$tip.label)))
	phylo$coalescentSurvivalProbability <- rep(1, (phylo$Nnode + length(phylo$tip.label)))
	class(phylo) <- c("binaryDatedTree", "phylo")
	return(phylo)
}


.extant.at.height <- function(h, tree)
{
	return( which( tree$heights <= h & tree$parentheights > h)  )
}


##########################
#CALCULATE INTERNAL STATES

.dPiAL <- function(h, y, parms, ...){
	# conditional on no coalescent
	p <- y[1:(parms$m-1)]
	A <- y[(parms$m):(2*parms$m-1)]
	L <- y[length(y)] #cumulative hazard
	pp <- c(p, 1 - sum(p))
	t <- parms$treeT - h
	
	if (USE_DISCRETE_FGY){
		if (FGY_INDEX < FGY_RESOLUTION &  h > FGY_H_BOUNDARIES[FGY_INDEX])  {FGY_INDEX <<- FGY_INDEX+1 ; return(.dPiAL(h, y, parms))}
		if ( FGY_INDEX > 1 ) { if ( h< FGY_H_BOUNDARIES[FGY_INDEX-1]) { FGY_INDEX <<- FGY_INDEX-1 ; return(.dPiAL(h, y, parms))} }
	}
	.G <- G.(t) 
	.F <- F.(t)
	.Y <- Y.(t) 
	X1 <- pmax(A / .Y, 0); X1[is.infinite(X1)] <- A[is.infinite(X1)]
	X2 <-  pmax( (.Y - A ) / .Y, 0); X2[is.infinite(X2)] <- A[is.infinite(X2)]
	
	FklXAk_Yk <- (.F * X1)
	if(FINITESIZECORRECTIONS){
		diag(FklXAk_Yk) <- diag(FklXAk_Yk) * .Y / (pmax(.Y,1.01)-1)
	}
	dL <- sum( FklXAk_Yk %*% X1 )
	
	#TODO would be faster to solve m- 1 equations since sum(A) is conserved
	dA <- c(+ as.vector(.G %*% X1) 
					- as.vector(colSums(.G) )* X1
					- as.vector( t(.F) %*% X2 ) * X1
					+ as.vector( .F %*% X1) * X2)
	
	R <- t(.G) / .Y + t(.F * X2) / .Y
	R <- R * (1- diag(length(.Y)))
	R <- R  -rowSums(R) * diag(length(.Y))
	R[is.nan(R)] <- 0
	
	dPi <-  (t(R) %*% pp)[1:(length(.Y)-1)]
	return(list( c(dPi, dA, dL) ))
}


.solve.Pi.and.AL <- function(h0, h1, A0, L0, tree)
{
	P0 				<- diag(tree$m)[,1:(tree$m-1)]
	if (tree$m <=2) 
		P0 			<- t(t(P0))
	parameters 		<- list(treeT = tree$maxSampleTime, m = tree$m)
	Pi1s 			<- c()
	for (i in 1:tree$m)
	{
		y0 			<- c( P0[i,], A0, L0)
		tryCatch({
					out0 <- ode(y=y0, times=c(h0, h1), func=.dPiAL, parms=parameters, method=INTEGRATIONMETHOD ) 
				}, error = function(e) browser() )
		Pi1 		<- abs(out0[nrow(out0),2:(1 + tree$m-1)])
		if( sum(Pi1) > 1) 
		{
			Pi1 	<- Pi1 / sum(Pi1)
		}
		Pi1s 		<- rbind(Pi1s, Pi1 )
	}
	Pi1s 			<- cbind(Pi1s, 1 - rowSums(Pi1s))
	A1 				<- out0[nrow(out0), (tree$m + 1):(2*tree$m-1 + 1)]
	L1 				<- out0[nrow(out0), ncol(out0)]
	return ( list(  unname(Pi1s), unname(A1), unname (L1) ) )
}

.calculate.internal.states <- function(tree, maxHeight=0){
	eventTimes <- unique( sort(tree$heights) )
	tree$maxHeight <- maxHeight
	if (maxHeight) { 
		eventTimes <- eventTimes[eventTimes<=maxHeight]}
	S <- 1
	L <- 0
	for (ih in 1:(length(eventTimes)-1)){
		h0 <- eventTimes[ih]
		h1 <- eventTimes[ih+1]
		#get A0, process new samples, calculate lstate for new lines
		extantLines <- .extant.at.height(h0, tree)
		if (length(extantLines) > 1 ){
			A0 <- colSums(tree$mstates[extantLines,])
		} else if (length(extantLines)==1){ A0 <- tree$mstates[extantLines,] }
		else{ browser() }
		
		#solve P[i,], nlft, and S=exp(-L)
		out <- .solve.Pi.and.AL(h0, h1, A0, L, tree)
		P <- out[[1]]
		A <- out[[2]]
		L <- out[[3]]
		
		#update midstates
		for (u in extantLines){
			tree$mstates[u,] <- t(P) %*% tree$mstates[u,]
			tree$mstates[u,] <- tree$mstates[u,] / sum(tree$mstates[u,])
		}
		
		#if applicable: update ustate & calculate lstate of new line 
		newNodes <- which( tree$heights == h1)
		newNodes <- newNodes[newNodes > length(tree$sampleTimes)]
		.F <- F.(tree$maxSampleTime - h1)
		.Y <- Y.(tree$maxSampleTime - h1)
		S <- exp(-L)
		
		#TODO option to return -Inf in this situation: 
		#if (sum(.Y) < length(extantLines) ) S <- 0 
		
		for (alpha in newNodes){
			u <- tree$daughters[alpha,1]
			v <- tree$daughters[alpha,2]
			tree$ustates[u,] <- tree$mstate[u,]
			tree$ustates[v,] <- tree$mstate[v,]
			
			FklXpuk_Yk <- (.F * tree$ustates[u,]/.Y)
			FklXpvk_Yk <- (.F * tree$ustates[v,]/.Y)
			# option for finite size corrections:
			if (FINITESIZECORRECTIONS){
				diag(FklXpuk_Yk) <- diag(FklXpuk_Yk) * .Y / (pmax(.Y,1.01)-1)
				diag(FklXpvk_Yk) <- diag(FklXpvk_Yk) * .Y / (pmax(.Y,1.01)-1)
				#TODO need additional finite size correction: update pvl conditional on state(u)==k
			}
			FklXpuk_Yk[is.nan(FklXpuk_Yk)] <- 0
			FklXpvk_Yk[is.nan(FklXpvk_Yk)] <- 0
			vk_Yk <- pmin(pmax(tree$ustates[v,]/.Y, 0),1); vk_Yk[is.nan(vk_Yk)] <- 0
			uk_Yk <- pmin(pmax(tree$ustates[u,]/.Y, 0),1); uk_Yk[is.nan(uk_Yk)] <- 0
			ratekl <- FklXpuk_Yk %*% vk_Yk + FklXpvk_Yk %*% uk_Yk
			
			tree$lstates[alpha,] <- ratekl / sum(ratekl)
			tree$mstates[alpha,] <- ratekl / sum(ratekl)
			tree$coalescentRates[alpha] <- sum(ratekl) 
			tree$coalescentSurvivalProbability[alpha] <- S
			
			# finite size corrections for lines not involved in coalescent
			if (FINITESIZECORRECTIONS)
			{
				p_a <- tree$lstates[alpha,]
				for (i in extantLines){
					if (i!=u & i!=v){
						p_i <- tree$mstates[i,]
						fterm <- p_i * p_a * pmax(.Y-1,0)/pmax(.Y-p_a, 1e-9)
						smat <-  t( matrix( rep( p_a * .Y/pmax(.Y-p_i,1e-9), tree$m), nrow=tree$m) )
						diag(smat) <- 0
						sterm <- p_i * rowSums(smat)
						tree$mstates[i,] <- fterm + sterm
						if (sum(is.na(fterm+sterm)) > 0) browser()
						if (sum(tree$mstates[i,]) <= 0) { tree$mstates[i,] <- .Y / sum(.Y) } 
						else{ tree$mstates[i,] <- tree$mstates[i,] / sum(tree$mstates[i,]) }
					}
				}
			}
		}
		
		if (length(newNodes) > 0) {L<-0}
	}
	return(tree)
}
#############################
.start.discrete.rates <- function(fgyResolution, maxHeight) 
{
	FGY_RESOLUTION		<<- fgyResolution
	#~ speed up calculation of FGY by discretizing & pre-caching
	F.bak 				<<- F.
	G.bak 				<<- G.
	Y.bak 				<<- Y.
	USE_DISCRETE_FGY 	<<- TRUE
	FGY_H_BOUNDARIES 	<<- seq(0, maxHeight, length.out = fgyResolution) 
	FGY_H_BOUNDARIES 	<<- FGY_H_BOUNDARIES + FGY_H_BOUNDARIES[2]/2
	FGY_INDEX 			<<- 1 #update in desolve:ode
	F_DISCRETE 			<<- lapply( FGY_H_BOUNDARIES, function(h) { F.(maxHeight-h) })
	G_DISCRETE 			<<- lapply( FGY_H_BOUNDARIES, function(h) { G.(maxHeight-h) })
	Y_DISCRETE 			<<- lapply( FGY_H_BOUNDARIES, function(h) { Y.(maxHeight-h) })
	F. 					<<- function(t) { F_DISCRETE[[FGY_INDEX]] } #note does not actually use arg t
	G. 					<<- function(t) { G_DISCRETE[[FGY_INDEX]] }
	Y. 					<<- function(t) { Y_DISCRETE[[FGY_INDEX]] }
}
#############################
.end.discrete.rates <- function(){
	#reset fgy functions
	F. <<- F.bak
	G. <<- G.bak
	Y. <<- Y.bak
}
#############################
# CALCULATE LIKELIHOOD
#' @export
coalescent.log.likelihood <- function(bdt, integrationMethod = 'rk4', finiteSizeCorrections=TRUE, maxHeight=0, discretizeRates=FALSE, fgyResolution = 100){
	print(paste(date(), 'start likelihood'))
	
	if (discretizeRates) { .start.discrete.rates(fgyResolution, maxHeight=max(bdt$heights) ) }
	else{ USE_DISCRETE_FGY <<- FALSE }
	
	# bdt : binaryDatedTree instance
	INTEGRATIONMETHOD <<- integrationMethod
	FINITESIZECORRECTIONS <<- finiteSizeCorrections
	tree <- .calculate.internal.states(bdt, maxHeight=maxHeight)
	i<- (length(tree$sampleTimes)+1):(tree$Nnode + length(tree$tip.label))
	if (maxHeight) { 
		internalHeights <- tree$heights[(length(tree$tip.label)+1):length(tree$heights)]
		i <- i[internalHeights <= maxHeight] }
	ll <- sum( log(tree$coalescentRates[i]) ) + sum( log(tree$coalescentSurvivalProbability[i]) )
	if (maxHeight) { ll<- tree$Nnode *  ll/length(i)}
	if (is.nan(ll) | is.na(ll) ) ll <- -Inf
	
	if (discretizeRates) { .end.discrete.rates()}
	print(paste(date(), 'finish likelihood', ll))
	return(ll)
}


#' Simulate binary dated tree
#' @export
simulatedBinaryDatedTree <- function( x, ...) UseMethod("simulatedBinaryDatedTree")
simulatedBinaryDatedTree.default <- function(sampleTime, sampleStates, discretizeRates=FALSE, fgyResolution = 100) 
{
	require(ape)
#~ simulates a coalescent tree, assumes F., G. and Y. are defined
#~ same attributes are defined as binaryDatedTree
#~ <preliminaries>
	n 			<- nrow(sampleStates)
	sampleTimes <- rep( sampleTime, n) 
	Nnode 		<-  n-1
	
	# NOTE when discretizing rates, this will neglect any changes in rates for t < 0
	if (discretizeRates) 
	{ 
		.start.discrete.rates(fgyResolution, maxHeight = max(sampleTimes)) 
	}
	else
	{ 
		USE_DISCRETE_FGY <<- FALSE 
	}
	#edge <- c()
	#edge.length <- c()
	edge.length 	<- rep(-1, Nnode + n-1) # should not have root edge
	edge			<- matrix(-1, (Nnode + n-1), 2)
	tip.label 		<- as.character( 1:n )
	maxSampleTime  	= T <- max(sampleTimes)
	heights 		<- rep(0, (Nnode + n) )
	parentheights 	<- rep(-1, (Nnode + n) )
	heights[1:n] 	<- maxSampleTime - sampleTimes
	inEdgeMap 		<- rep(-1, Nnode + n)
	outEdgeMap 		<- matrix(-1, (Nnode + n), 2)
	parent 			<- 1:(Nnode + n) 
	daughters 		<- matrix(-1, (Nnode + n), 2)
	m	 			<- dim(sampleStates)[2]
	lstates 		<- matrix(-1, (Nnode + n), m)
	mstates 		<- matrix(-1, (Nnode + n), m)
	ustates 		<- matrix(-1, (Nnode + n), m)
	lstates[1:n,] 	<- sampleStates
	mstates[1:n,] 	<- sampleStates
	ustates[1:n,] 	<- sampleStates
#~ </preliminaries>
	
#~ <survival time to next event>
	calculate.rates <- function(h, parms)
	{
		# eliminate diag elements for migration
		t 					<- parms$maxSampleTime - h
		.G 					<- G.(t) 
		.F 					<- F.(t)
		.Y 					<- Y.(t) 
		X1 					<- pmax( parms$A / .Y, 0)
		X2 					<-  pmax( (.Y - parms$A ) / .Y, 0)
		X1[is.nan(X1)] 		<- 0
		X2[is.nan(X2)] 		<- 0
		X1[is.infinite(X1)] <- 0
		X2[is.infinite(X2)] <- 0
		X1mat 				<- matrix(X1,nrow=parms$m,ncol=parms$m, byrow=TRUE)
		tX1mat 				<- matrix(X1,nrow=parms$m,ncol=parms$m, byrow=FALSE)
		tX2mat 				<- matrix(X2, nrow=parms$m,ncol=parms$m, byrow=FALSE)
		lambdaCoalescent 	<-  .F * X1mat * tX1mat
		lambdaMigration 	<-  t(.G * X1mat  )
		lambdaInvisibleTransmission	<-  t(.F * tX2mat * X1mat)
		diag(lambdaMigration) 		<- 0
		diag(lambdaInvisibleTransmission) <- 0
		#browser()
		return(list( X1=X1, X2=X2, lambdaCoalescent=lambdaCoalescent, lambdaMigration=lambdaMigration, lambdaInvisibleTransmission=lambdaInvisibleTransmission) )
	}
	dtheta <- function(h, theta, parms, ...){
		if (is.nan(theta)) return(list(NaN))
		if (theta <= parms$r) return(list(NaN)) #terminate early if survival time found
		rates <- calculate.rates(h, parms)
		return( list( -theta * (sum(rates$lambdaCoalescent) + sum(rates$lambdaMigration) + sum(rates$lambdaInvisibleTransmission)) ) )
	}
#~ </survival time to next event>
	
#~ <simulate tree>
#TODO first version will work for homochronous, then adapt for account for mult samp times
#TODO these trees are still biased; check rate matrices
	extant <- 1:n
	lineageCounter <- n+1 # next lineage will have this index
	A <- colSums( mstates[extant,])
	h <- 0
	notdone <- TRUE
	lastExtant <- n #DEBUG
	while(notdone){
		r <- runif(1)
		theta0 <- 1
		parms <- list(A=A, m=m, maxSampleTime=maxSampleTime, r=r)
		
		# find appropriate time axis
		rates <- calculate.rates( h, parms) 
		parms$rates <- rates
		lambda <- (sum(rates$lambdaCoalescent) + sum(rates$lambdaMigration) + sum(rates$lambdaInvisibleTransmission))
		
		# solve survival time
		# NOTE more accurate to interpolate approxfun( thetaTimes, o[,2] ), & invert to find o[o[,2]==r,1]
		if (discretizeRates){
			eventTime <- rexp(1, rate=lambda) + h
			while (FGY_INDEX < FGY_RESOLUTION &  eventTime > FGY_H_BOUNDARIES[FGY_INDEX]) {
				if (FGY_INDEX < FGY_RESOLUTION &  eventTime > FGY_H_BOUNDARIES[FGY_INDEX])  {
					eventTime <- FGY_H_BOUNDARIES[FGY_INDEX] 
					FGY_INDEX <<- FGY_INDEX+1 ; 
					tryCatch({
								rates <- calculate.rates( eventTime, parms) ; parms$rates <- rates
							}, error = function(e) browser() )
					lambda <- (sum(rates$lambdaCoalescent) + sum(rates$lambdaMigration) + sum(rates$lambdaInvisibleTransmission))
					eventTime <- eventTime + rexp(1, rate=lambda)
				}
			}
		} else{
			hub <- h - log(1e-6)/lambda
			thetaTimes <- seq(h, hub, length.out = SIMULATIONTIMERESOLUTION)
			tryCatch( { o <- ode(y=c(theta0), times=thetaTimes, func=dtheta, parms = parms)}, error =  function(e) {browser() } )
			eventTimeIndex <- match( NaN, o[,2])-1
			if (is.na(eventTimeIndex) ) eventTimes <- thetaTimes[length(thetaTimes)] # NOTE more accurate to continue integration in this case
			eventTime <- thetaTimes[eventTimeIndex]
		}
		
		# which event? 2m2 possibilities
		cumulativeRatesVector <- cumsum( c( as.vector(rates$lambdaCoalescent), as.vector(rates$lambdaMigration+rates$lambdaInvisibleTransmission)) )
		# as.vector unravels matrices by column
		cumulativeRatesVector <- cumulativeRatesVector / cumulativeRatesVector[length(cumulativeRatesVector) ]
		rwhichevent <- runif(1)
		indexEvent <- match(TRUE, cumulativeRatesVector > rwhichevent)
		if (is.na(indexEvent)){ 
			warning('Rates at ', eventTime, 'are NA or all rates are zero')
			indexEvent <- round( runif(1) * length(cumulativeRatesVector) )
		}
		cumulativeRatesVector
		indexEvent
		
		if (indexEvent <= m*m){ #coalescent event
			l <- 1+floor( (indexEvent-1) / m )
			k <-  1+ ( (indexEvent-1) %% m) #1+ ((indexEvent+1) %% m)
			kextant <- extant[ mstates[extant, k]==1 ]
			lextant <- extant[ mstates[extant, l]==1 ]
			if (k==l){ # guard against A[k]<2
				if (length(kextant) >=2) {ij <- sample(kextant, 2) }
				else { ij <- sample(extant, 2) }
			} else{
				# horribly ugly code thanks to stupid default behavior of 'sample':
				if ((length(kextant)>1) & (length(lextant)>1)) { ij <- c( sample(kextant, 1), sample(lextant, 1) ) }
				else if (length(kextant)==0 | length(lextant)==0) { warning('kextant length zero'); ij <- sample(extant, 2) }
				else if (length(kextant)==1 & length(lextant)==1) { ij <- c(kextant, lextant) }
				else if (length(kextant)==1) { ij <- c(kextant, sample(lextant, 1)) }
				else if (length(lextant)==1) { ij <- c(sample(kextant, 1), lextant) }
				else { browser() }
			}
			# set new lineage
			extant <- c(extant, lineageCounter)
			extant <- extant[extant!=ij[1] & extant!=ij[2]]
			heights[lineageCounter] <- eventTime
			lstates[lineageCounter,] <- diag(m)[k,]
			mstates[lineageCounter,] <- diag(m)[k,]
			inEdgeMap[ij] <- lineageCounter
			outEdgeMap[lineageCounter,] <- ij
			parent[ij] <- lineageCounter
			parentheights[ij] <- eventTime
			daughters[lineageCounter,] <- ij
			#edge <- rbind(edge, c(lineageCounter, ij[1]))
			#edge.length <- c(edge.length, eventTime - heights[ij[1]] )
			#edge <- rbind(edge, c(lineageCounter, ij[2]))
			#edge.length <- c(edge.length, eventTime - heights[ij[2]] )
			edge[ij[1],] <- c(lineageCounter, ij[1]) # 
			edge.length[ij[1]] <- eventTime - heights[ij[1]] 
			edge[ij[2],] <- c(lineageCounter, ij[2]) # edge <- rbind(edge, c(lineageCounter, ij[2]))
			edge.length[ij[2]] <- eventTime - heights[ij[2]] 
			lineageCounter <- lineageCounter+1
			#print(paste(date(), 'coalescent', h, eventTime,  k, l, length(extant)))
			#if (length(extant)>= lastExtant) browser() #DEBUG
			#lastExtant <- length(extant) #DEBUG
			#if (rates$lambdaCoalescent[k,l]==0) {print('coalescent'); browser()}
		} else{ #migration event
			indexEvent <- indexEvent- m*m
			l <- 1+floor( (indexEvent-1) / m )
			k <- 1+ ( (indexEvent-1) %% m) #1+ ( (indexEvent+1) %% m)
			kextant <- extant[ mstates[extant, k]==1 ]
			if (length(kextant) > 1) { i <- sample(kextant, 1) }
			else { i <- sample(extant, 1) }
			mstates[i,] <- diag(m)[l,]
			#print(paste( date(),'migration', h, eventTime, k, l) )
			#print(paste(  k, l) )
			#migrations <- rates$lambdaMigration+rates$lambdaInvisibleTransmission
			#if (migrations[k,l]==0) {print('migrations'); browser()}
		}
		if(length(extant) <= 1 ) {notdone <- FALSE}
		else{
			A <- colSums( mstates[extant,])
			h <-  eventTime
		}
	}
#~ </simulate tree>
	if (discretizeRates) { .end.discrete.rates()}
#~ <assemble class>
	self<- list(edge=edge, edge.length=edge.length, Nnode=Nnode, tip.label=tip.label, heights=heights, parentheights=parentheights, parent=parent, daughters=daughters, lstates=lstates, mstates=mstates, ustates=ustates, m = m, sampleTimes = sampleTimes, sampleStates= sampleStates, maxSampleTime=maxSampleTime, inEdgeMap = inEdgeMap, outEdgeMap=outEdgeMap)
	class(self) <- c("binaryDatedTree", "phylo")
#~ </>
	
	
#~ <reorder edges for compatibility with ape::phylo functions> 
#~ (ideally ape would not care about the edge order, but actually most functions assume a certain order)
	sampleTimes2 <- sampleTimes; names(sampleTimes2) <- tip.label
	sampleStates2 <- sampleStates; rownames(sampleStates2) <- tip.label
	phylo <- read.tree(text=write.tree(self) )
	sampleTimes2 <- sampleTimes2[phylo$tip.label]; sampleTimes2 <- unname(sampleTimes2)
	sampleStates2 <- sampleStates2[phylo$tip.label,]; sampleStates2 <- unname(sampleStates2)
	binaryDatedTree(phylo, sampleTimes2, sampleStates2)
#~ </>
}

#~ setOldClass("binaryDatedTree")
#~ plot.binaryDatedTree <- function(x,...) { 
#~ 	plot.phylo( read.tree(text=write.tree(x))  , ...)
#~ }
#~ setMethod("plot.phylo", signature(x="binaryDatedTree"), plot.binaryDatedTree)



#######################################
# SUMMARY STATISTIC GENERATOR
# TODO adapt for heterochronous samples
#~ m + m2+ m3
#' Calculate cluster size moments from model
#' @param sampleTime 
#' @param sampleStates	
#' @param maxTime
#' @param minTime		
#' @param timeResolution
#' @param discretizeRates
#' @param fgyResolution
#' @param integrationMethod
#' @export
#'
calculate.cluster.size.moments.from.model <- function(sampleTime, sampleStates , maxTime=NA, minTime = NA, timeResolution = 50, discretizeRates=FALSE, fgyResolution = 100 , integrationMethod = 'adams')
{
	require(deSolve)
	# assumes that F, G, Y are defined
	n 		<- nrow(sampleStates) 
	m 		<- ncol(sampleStates)
	sampleTimes <- rep(sampleTime, n)
	if (is.na(maxTime)) maxTime <- sampleTime 
	if (is.na(minTime)) minTime <- 0
	heights <- seq( 0, maxTime-minTime, length.out = timeResolution)
	

	if (discretizeRates) 
	{ 
		.start.discrete.rates(fgyResolution, maxHeight = max(sampleTimes)) 
	}
	else
	{ 
		USE_DISCRETE_FGY <<- FALSE 
	}
	INTEGRATIONMETHOD <<- integrationMethod
	
	
#~ 	#solve first aggregated moment X1
#~ 	# X_l^j : number type l descended from type j lineages
#~ 	X1 <- list() # index by tip type, value = heights X ancestorType matrix
#~ 	X1_interps <- list() #list of lists of approxfun, index by [[tip type]][[ancestor type]]
#~ 	dX1l <- function(h, X1l, parms, ...){ #same derivative function for all tip types
#~ 		t <- parms$maxTime - h
#~ 		.Y <- Y.(t)
#~ 		.F <- F.(t)
#~ 		.G <- G.(t)
#~ 		X1l_Y <- X1l / .Y
#~ 		.F %*% X1l_Y + .G %*% X1l_Y - colSums( .F + .G ) * X1l_Y
#~ 	}
	
	dXiXjXij <- function(h, XiXjXij, parms, ...)
	{
		#X are vectors indexed by ancestor type
		Xi 	<- XiXjXij[1:m]
		Xj 	<- XiXjXij[(m+1):(2*m)]
		Xij <- XiXjXij[(2*m+1):(length(XiXjXij))]
		t 	<- parms$maxTime - h
		if (USE_DISCRETE_FGY)
		{
			if (FGY_INDEX < FGY_RESOLUTION &  h > FGY_H_BOUNDARIES[FGY_INDEX])  
			{
				FGY_INDEX <<- FGY_INDEX+1 
				return(dXiXjXij(h, XiXjXij, parms))
			}
			if ( FGY_INDEX > 1 ) 
			{ 
				if ( h< FGY_H_BOUNDARIES[FGY_INDEX-1]) 
				{ 
					FGY_INDEX <<- FGY_INDEX-1 
					return(dXiXjXij(h, XiXjXij, parms))
				} 
			}
		}
		.Y 		<- Y.(t)
		.F 		<- F.(t)
		.G 		<- G.(t)
		Xi_Y 	<- Xi / .Y
		Xj_Y 	<- Xj / .Y
		Xij_Y	<- Xij / .Y
		csFpG 	<- colSums( .F + .G )
		dXi 	<- .F %*% Xi_Y + .G %*% Xi_Y - csFpG * Xi_Y
		dXj 	<- .F %*% Xj_Y + .G %*% Xj_Y - csFpG * Xj_Y
		dXij	<- .F %*% Xij_Y + .G %*% Xij_Y - csFpG * Xij_Y  + (.F %*% Xj_Y )*Xi_Y  + (.F %*% Xi_Y ) * Xj_Y
		list( c( dXi, dXj, dXij ) )
	}
	
	dA <- function(h, A, parms, ...)
	{
		if (USE_DISCRETE_FGY)
		{
			if (FGY_INDEX < FGY_RESOLUTION &  h > FGY_H_BOUNDARIES[FGY_INDEX])  
			{
				FGY_INDEX <<- FGY_INDEX+1 
				return(dA(h, A, parms))
			}
			if ( FGY_INDEX > 1 ) 
			{ 
				if ( h< FGY_H_BOUNDARIES[FGY_INDEX-1]) 
				{ 
					FGY_INDEX <<- FGY_INDEX-1 
					return(dA(h, A, parms))
				} 
			}
		}
		t 		<- parms$maxTime - h
		.Y 		<- Y.(t)
		.F 		<- F.(t)
		.G 		<- G.(t)
		A_Y 	<- A / .Y
		csFpG	<- colSums( .F + .G )
		list( .G %*% A_Y - csFpG * A_Y + (.F %*% A_Y) * pmax(1-A_Y, 0) )
	}
	
#~ 	for (l in 1:m){ #type of tip
#~ 	X1l_h0 <- diag(m)[l,] * A0 #vector corresponding to type of ancestor
#~ 	parameters <- list(maxTime = maxTime) 
#~ 	o<- ode(y=X1l_h0, times=heights, func=dX1l, parms=parameters, method=INTEGRATIONMETHOD ) 
#~ 		X1lk_interps <- list()
#~ 		for (k in 1:m){ 
#~ 		X1lk_interps[[k]] <- approxfun( o[,1], o[,(1+k)], yleft=o[1,(1+k)], yright=o[nrow(o),(1+k)] )
#~ 		}
#~ 	X1[[l]] <- o[,2:ncol(o)]
#~ 	X1_interps[[l]] <- X1lk_interps
#~ 	}
	
	# solve for A
	A0 			<- colSums(sampleStates) 
	parameters 	<- list(maxTime = maxTime)
	o 			<- ode(y=A0, times=heights, func=dA, parms=parameters, method=INTEGRATIONMETHOD ) 
	FGY_INDEX 	<<- 1
	A 			<- o[,2:(ncol(o))]
	
	# solve Xi Xj and Xij at same time, 3m variables, avoids approxfun nastiness
	# solve second aggregated moment X2 & moments
	Mij_h 		<- array(0, dim = c(m, m, timeResolution) )
	Mi_h		<- matrix(0, m, timeResolution)
	for (i in 1:m)
	{ # first tip type
		for (j in i:m)
		{ #second tip type
			#if (i!=j){
			parameters 		<- list(maxTime = maxTime) 
			Xij_h0 			<- rep(0, m)
			if(i==j) 
				Xij_h0[i]	<- A0[i]
			XiXjXij_h0 		<- c( diag(m)[i,] * A0, diag(m)[j,] * A0, Xij_h0) # each X is vector m(ancestor type) X 1
			o 				<- ode(y=XiXjXij_h0, times=heights, func=dXiXjXij, parms=parameters, method=INTEGRATIONMETHOD )
			# moments
			Xij_h 			<- o[, (1 + 2*m+1) : ncol(o) ] 
			if(i==j)
			{
				Xi_h		<- o[, 1+1:m]	
				tryCatch( 
					Mi_h[i,]<- rowSums( Xi_h   ) / rowSums(A)	# why avg over cols? / total lineages								
							, error = function(e) browser() )
			}		
			# covariances
			tryCatch( 					
				Mij_h[i,j,]	<- rowSums( Xij_h   ) / rowSums(A) 	# timeResolution X 1
				, error = function(e) browser() )
			Mij_h[j,i,] 	<- Mij_h[i,j,]			
		}
	}
	
	if (discretizeRates) 
		.end.discrete.rates()
	list( heights=heights, M= Mij_h, E=Mi_h, A = A)
}


#' Calculate cluster size moments from tree
#' @param bdt		binary dated tree
#' @param heights	vector numeric, heights at which to calculate cluster sizes	
#' @export
#'
calculate.cluster.size.moments.from.tree <- function(bdt, heights)
{
# bdt : binaryDatedTree
# heights : vector numeric, heights at which to calculate cluster sizes
	m 		<- ncol( bdt$sampleStates )
	Mi		<- matrix(0, m, length(heights))
	Mij 	<- array(0, dim=c(m, m, length(heights)) )
	is.tip 	<- function(u) 
	{ 
		ifelse(is.na(bdt$daughters[u,1]), TRUE, FALSE)
	}
	identify.tips <- function(u) 
	{
		if (is.tip(u)) return(u)
		uv	<- bdt$daughters[u,]
		c(identify.tips(uv[1]), identify.tips(uv[2]) ) 
	}
	calculate.moment2 <- function(u, i, j) 
	{
		tips	<- identify.tips(u)
		itips 	<-  sum( bdt$sampleStates[tips,i] )
		jtips 	<-  sum( bdt$sampleStates[tips,j] )
		itips * jtips
	}
	calculate.mean.moment2 <- function(extant, i, j)
	{
		if (length(extant)==0) return(NA)
		mean( sapply( extant, function(u) calculate.moment2(u, i, j) ) ) 
	}
	for (ih in 1:length(heights))
	{
		h 		<- heights[ih]
		extant 	<- .extant.at.height(h, bdt) 
		if (length(extant) <=1)
		{
			for (i in 1:m)
			{
				for (j in i:m)
				{
					Mij[i,j,ih] <- Mij[i,j,ih-1]
					Mij[j,i,ih] <- Mij[i,j,ih]
				}
			}
		} 
		else
		{
			for (i in 1:m)
			{
				for (j in i:m)
				{
					Mij[i,j,ih] <- calculate.mean.moment2(extant,i,j)
					Mij[j,i,ih] <- Mij[i,j,ih]
				}
			}
		}
	}	
	Mij
}
