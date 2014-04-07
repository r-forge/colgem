#' this file contains all R functions of the rcolgem package
#' @import ape
#' @import deSolve

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



.calculate.internal.states.unstructuredModel <- function(tree, maxHeight=0, globalParms = NULL)
{
	if (!is.null(globalParms)){
		F. <- globalParms$F.; G. <- globalParms$G.; Y. <- globalParms$Y. ; 
	   USE_DISCRETE_FGY <- globalParms$USE_DISCRETE_FGY; 
	   FGY_INDEX <- globalParms$FGY_INDEX
	   FGY_RESOLUTION <- globalParms$FGY_RESOLUTION
	   FGY_H_BOUNDARIES <- globalParms$FGY_H_BOUNDARIES
	   FINITESIZECORRECTIONS <- globalParms$FINITESIZECORRECTIONS
	}
	eventTimes <- unique( sort(tree$heights) )
	tree$maxHeight <- maxHeight
	if (maxHeight) { 
		eventTimes <- eventTimes[eventTimes<=maxHeight]}
	S <- 1
	L <- 0
	
	FGY_INDEX <- 1 
	
	get.fgy <- function(h,t=NULL)
	{
		if (globalParms$USE_DISCRETE_FGY)
		{
			while (FGY_INDEX < globalParms$FGY_RESOLUTION &  h > globalParms$FGY_H_BOUNDARIES[FGY_INDEX])  
			{
				FGY_INDEX <- FGY_INDEX+1 ; 
			}
			if ( (FGY_INDEX > 1) ) 
			{
				while ( h< globalParms$FGY_H_BOUNDARIES[FGY_INDEX-1]) 
				{
					if (FGY_INDEX==1) break;
					FGY_INDEX <- FGY_INDEX-1 ; 
					#return(dA(h, A, parms))
				} 
			}
			.Y		<- globalParms$Y.(FGY_INDEX)
			.F		<- globalParms$F.(FGY_INDEX)
			.G		<- globalParms$G.(FGY_INDEX)
		} else
		{
			.Y		<- globalParms$Y.(t)
			.F		<- globalParms$F.(t)
			.G		<- globalParms$G.(t)
			FGY_INDEX <- NA
		}
		list(.F=.F, .G=.G, .Y=.Y, FGY_INDEX=FGY_INDEX)
	}
	
	#DEBUG
	#with( get.fgy(0), {print(.Y); print(.F)  })
	#browser()
	
	.dL.unstructuredModel <- function(h, L, parms, FGY_INDEX, ...){
		# conditional on no coalescent
		t <- parms$treeT - h
		A <- parms$A
		fgy <- get.fgy(h,t)
		with( fgy, 
		{
			.F * (A / .Y)^2# *  (A-1) / max((.Y-1), 0.01)
			#.F * (A / .Y) *  (A-1) / max((.Y-1), 0.01)
		}) -> dL 
		return( list(dL, FGY_INDEX=fgy$FGY_INDEX) )
	}
	.solve.L <- function(h0, h1, L,  A0 ) 
	{
		parameters 		<- list(treeT = tree$maxSampleTime, m = tree$m, A=A0)
			tryCatch({
						out0 <- ode(y=L, times=c(h0, h1), func=.dL.unstructuredModel, parms=parameters, FGY_INDEX=FGY_INDEX, method=globalParms$INTEGRATIONMETHOD ) 
					}, error = function(e) browser() )
			L1 <- unname( out0[nrow(out0),2] )
			FGY_INDEX <- unname( out0[nrow(out0),3] )
			return( list( L1, FGY_INDEX) )
	}
	
	
	
	for (ih in 1:(length(eventTimes)-1)){
		h0 <- eventTimes[ih]
		h1 <- eventTimes[ih+1]
		#get A0, process new samples, calculate lstate for new lines
		extantLines <- .extant.at.height(h0, tree)
		A0 <- length(extantLines)
		
		o <- .solve.L(h0, h1, L,  A0 ) 
		L <- o[[1]]
		FGY_INDEX <- o[[2]]
		
		newNodes <- which( tree$heights == h1)
		newNodes <- newNodes[newNodes > length(tree$sampleTimes)]
		if (length(newNodes) > 0)
		{
			fgy<- get.fgy(h1,  bdt$maxSampleTime - h1 )
			.Y <- fgy$.Y
			.F <- fgy$.F
			FGY_INDEX <- fgy$FGY_INDEX
			
			S <- exp(-L)
			
			#TODO option to return -Inf in this situation: 
			#tryCatch( { if (sum(.Y) < length(extantLines) ) {S <- 0; } }, error = function(e) browser() )
			#print(paste(sum(.Y),length(extantLines)));
			
			for (alpha in newNodes){
				X1 <- A0 / .Y
				# option for finite size corrections:
				#if (FINITESIZECORRECTIONS){
				#	X1 <- X1 * .Y / (max(.Y,1.01)-1)
				#}
				
				tree$lstates[alpha,] <- 1
				tree$mstates[alpha,] <- 1
				ratekl <- .F * (A0 / .Y)^2# *  (A0-1) / max((.Y-1), 0.01)
				#ratekl <- .F * (A0 / .Y) *  (A0-1) / max((.Y-1), 0.01)
				tree$coalescentRates[alpha] <- ratekl
				tree$coalescentSurvivalProbability[alpha] <- S
			}
			
			L<-0
		}
	}
	return( tree )
}


.calculate.internal.states <- function(tree, maxHeight=0, globalParms=NULL){
print('testing version')
	if (!is.null(globalParms)){
		F. <- globalParms$F.; G. <- globalParms$G.; Y. <- globalParms$Y. ; 
	   USE_DISCRETE_FGY <- globalParms$USE_DISCRETE_FGY; 
	   FGY_INDEX <- globalParms$FGY_INDEX
	   FGY_RESOLUTION <- globalParms$FGY_RESOLUTION
	   FGY_H_BOUNDARIES <- globalParms$FGY_H_BOUNDARIES
	}
	
	FGY_INDEX <- 1
	get.fgy <- function(h,t=NULL)
	{
		if (globalParms$USE_DISCRETE_FGY)
		{
			FGY_INDEX <- min( 1+floor( globalParms$FGY_RESOLUTION * h / globalParms$maxHeight ), globalParms$FGY_RESOLUTION)
			.Y		<- globalParms$Y.(FGY_INDEX)
			.F		<- globalParms$F.(FGY_INDEX)
			.G		<- globalParms$G.(FGY_INDEX)
			
		} else
		{
			.Y		<- globalParms$Y.(t)
			.F		<- globalParms$F.(t)
			.G		<- globalParms$G.(t)
		}
		list(.F=.F, .G=.G, .Y=.Y,  FGY_INDEX=FGY_INDEX)
	}
	
	
	.dQ.L6 <- function(h, y, parms, globalParms=NULL, ...){
		# see notes march 25
		Q <- pmax(matrix(y[1:parms$m^2], nrow=parms$m, ncol=parms$m), 0) #
		.omc <- colSums( Q ) # 'one minus c'
		P <- t( t(Q) / .omc ) #col wise
		L <- y[length(y)] #cumulative hazard
		t <- parms$treeT - h
		
		A <- as.vector( P %*% parms$A0 )
		
		with(get.fgy(h, t), {
			X1 <- pmax(A / .Y, 0); X1[is.infinite(X1)] <- 1 #A[is.infinite(X1)]
			X2 <- pmax( ((A-1)/(.Y-1)), 0);  X2[is.infinite(X2)] <- 1
						
			F1 <- .F; diag(F1) <- 0
			F2 <- diag( diag(.F) )
			dL <- X1 %*% F1 %*% X1 
			    + X1 %*% F2 %*% X2
			
			# non-coalescent transitions:
			O <- t( .G  + .F ) / .Y 
			diag(O) <- 0
			diag(O) <- -rowSums(O)
			
			# loss due to coalescent
			C <- diag( colSums( X1*(.F + t(.F)) ) / .Y )  #transmissions by index resulting in coalescent
			
			# final
			R <- O-C
			R[is.nan(R)] <- 0
			
			dQ <- t(R) %*% P 
			dQ[is.nan(dQ)] <- 0
			return(list( c(dQ, as.vector( dL) ) ))
		})
	}
	
	.solve.P.L3 <- function(h0, h1, A0, L0, tree, globalParms=NULL)
	{
		Q0 <- diag(tree$m)
		parameters 		<- list(treeT = tree$maxSampleTime, m = tree$m, A0 = A0)
		y0 <- c( as.vector( Q0), L0 ) #TODO
		o <- ode(y=y0, times=c(h0, h1), func=.dQ.L6, parms=parameters, globalParms=globalParms, method=globalParms$INTEGRATIONMETHOD ) 
		Q1 		<- t( matrix(  abs(o[nrow(o),2:(1 + tree$m^2)]) , nrow=tree$m) ) #NOTE the transpose
		P1 <- Q1 / rowSums(Q1) # NOTE renormalization
		if (sum(is.nan(P1))>0) browser()
		A1 <- t(P1) %*% A0
		L1 <- o[nrow(o), ncol(o)]
		return ( list(  unname(P1), unname(A1), unname (L1) ) )
	}
	
	
	
	eventTimes <- unique( sort(tree$heights) )
	tree$maxHeight <- maxHeight
	if (maxHeight) { 
		eventTimes <- eventTimes[eventTimes<=maxHeight]}
	if (max(eventTimes) > tree$maxSampleTime)  warning('Node heights of tree extend to negative time axis. Results may not be accurate. Try rescaling time axis and sample times. ')
	S <- 1
	L <- 0
	hAs <- c() #debug
	for (ih in 1:(length(eventTimes)-1)){
		h0 <- eventTimes[ih]
		h1 <- eventTimes[ih+1]
		#get A0, process new samples, calculate lstate for new lines
		extantLines <- .extant.at.height(h0, tree)
		if (length(extantLines) > 1 ){
			A0 <- colSums(tree$mstates[extantLines,])
		} else if (length(extantLines)==1){ A0 <- tree$mstates[extantLines,] }
		else{ browser() }
		hAs <- rbind( hAs, c(h0, A0)) #debug
		
		out <- .solve.P.L3(h0, h1, A0, L, tree, globalParms=globalParms)
		P <- out[[1]]
		A <- out[[2]]
		L <- out[[3]]
		
		#update midstates
		tree$mstates[extantLines,] <- t( t(P) %*% t(tree$mstates[extantLines,])  )
		tree$mstates[extantLines,] <- tree$mstates[extantLines,] / rowSums(tree$mstates[extantLines,])
		
		#if applicable: update ustate & calculate lstate of new line 
		newNodes <- which( tree$heights == h1)
		newNodes <- newNodes[newNodes > length(tree$sampleTimes)]
		
		fgy <- get.fgy(h1, globalParms$maxHeight-h1)
		with(fgy, 
		{
			S <- exp(-L)
			
			#TODO option to return -Inf in this situation: 
			if (sum(.Y) < length(extantLines) ) S <- 0 
			
			for (alpha in newNodes){
				u <- tree$daughters[alpha,1]
				v <- tree$daughters[alpha,2]
				
				{
					tree$ustates[u,] <- tree$mstate[u,]
					tree$ustates[v,] <- tree$mstate[v,]
					
					FklXpuk_Yk <- (.F * tree$ustates[u,]/.Y)
					FklXpvk_Yk <- (.F * tree$ustates[v,]/.Y)
					FklXpuk_Yk[is.nan(FklXpuk_Yk)] <- 0
					FklXpvk_Yk[is.nan(FklXpvk_Yk)] <- 0
					vk_Yk <- pmin(pmax(tree$ustates[v,]/.Y, 0),1); vk_Yk[is.nan(vk_Yk)] <- 0
					uk_Yk <- pmin(pmax(tree$ustates[u,]/.Y, 0),1); uk_Yk[is.nan(uk_Yk)] <- 0
					ratekl <- FklXpuk_Yk %*% vk_Yk + FklXpvk_Yk %*% uk_Yk
					
					rate_uv <- uk_Yk %*% .F %*% vk_Yk
					rate_vu <- vk_Yk %*% .F %*% uk_Yk
					putrans <- rate_uv / (rate_uv + rate_vu)
					if (is.nan(putrans)) putrans <- .5
					pvtrans <- 1 - putrans
				}
				
				if (sum(ratekl)==0) ratekl <- rep(1/tree$m, tree$m) * 1e-6 # so that line states are not nan
				
				# ancestor state at new node
					tree$lstates[alpha,] <- ratekl / sum(ratekl)
					tree$mstates[alpha,] <- ratekl / sum(ratekl)
				# for likelihood:
				tree$coalescentRates[alpha] <- sum(ratekl) 
				tree$coalescentSurvivalProbability[alpha] <- S
				tree$logCoalescentSurvivalProbability[alpha] <- -L
				
				
				# finite size corrections for lines not involved in coalescent
				# can be important when sample fraction is large and/or population size is small
				{ 
				#highly vectorized version- should be fast
				# modification from Genetics '12: uses A instead of Y
					p_a     <- tree$lstates[alpha,]
					p_i_mat <- tree$mstates[extantLines,]
					A_mat   <- t( matrix( A, nrow=m, ncol=length(extantLines) )  )
					p_a_mat <- t( matrix(p_a, nrow=m, ncol=length(extantLines)) )
					rho_mat <- A_mat /  (A_mat - p_i_mat)
					rterms  <- p_a_mat / (A_mat - p_i_mat)
					lterms  <- rho_mat %*% p_a
					lterms <- matrix(lterms, nrow=length(lterms), ncol=m) #maybe this is faster than sapply
					new_p_i <- p_i_mat * (lterms - rterms)
					new_p_i <- new_p_i / rowSums(new_p_i)
					tree$mstates[extantLines,] <- new_p_i
				}
			}
			tree
		}) -> tree
		if (length(newNodes) > 0) {L<-0}
	}
#~ 	return(tree)
return( list(tree, hAs)) #debug
} 


.start.discrete.rates <- function(fgyResolution, maxSampleTime, globalParms) 
{ #version that does not generate global variables
	F. <- globalParms$F.; G. <- globalParms$G.; Y. <- globalParms$Y. ; 
	globalParms$FGY_RESOLUTION		<- fgyResolution
	#~ speed up calculation of FGY by discretizing & pre-caching
	globalParms$F.bak 				<- F.
	globalParms$G.bak 				<- G.
	globalParms$Y.bak 				<- Y.
	globalParms$USE_DISCRETE_FGY 	<- TRUE
	globalParms$maxHeight			<- maxSampleTime #NOTE: This does not allow FGY defined for t<0
	globalParms$FGY_H_BOUNDARIES 	<- seq(0, maxSampleTime, length.out = fgyResolution) 
#~ 	globalParms$FGY_H_BOUNDARIES 	<- globalParms$FGY_H_BOUNDARIES + globalParms$FGY_H_BOUNDARIES[2]/2 #TODO is there a better way to do this offset? 
	globalParms$FGY_INDEX 			<- 1 #update in desolve:ode
	globalParms$F_DISCRETE 			<- lapply( globalParms$FGY_H_BOUNDARIES, function(h) { F.(maxSampleTime-h) })
	globalParms$G_DISCRETE 			<- lapply( globalParms$FGY_H_BOUNDARIES, function(h) { G.(maxSampleTime-h) })
	globalParms$Y_DISCRETE 			<- lapply( globalParms$FGY_H_BOUNDARIES, function(h) { Y.(maxSampleTime-h) })
	globalParms$F. 					<- function(FGY_INDEX) { globalParms$F_DISCRETE[[FGY_INDEX]] } #note does not actually use arg t
	globalParms$G. 					<- function(FGY_INDEX) { globalParms$G_DISCRETE[[FGY_INDEX]] }
	globalParms$Y. 					<- function(FGY_INDEX) { globalParms$Y_DISCRETE[[FGY_INDEX]] }
	globalParms
}

#############################
.end.discrete.rates <- function(globalParms){
	#reset fgy functions
	# does not use global variables
	globalParms$F. <- globalParms$F.bak
	globalParms$G. <- globalParms$G.bak
	globalParms$Y. <- globalParms$Y.bak
	globalParms
}
#############################

# CALCULATE LIKELIHOOD
#' @export
coalescent.log.likelihood <- function(bdt, FGY=NULL, integrationMethod = 'rk4',  maxHeight=0, discretizeRates=TRUE, fgyResolution = 200)
{
	#print(paste(date(), 'start likelihood'))
	globalParms <- list( USE_DISCRETE_FGY = discretizeRates , INTEGRATIONMETHOD=integrationMethod )
	if (!is.null(FGY)) { globalParms <- modifyList( globalParms, FGY) }
	else {globalParms <- modifyList( globalParms, list( F. = F., G. = G., Y. = Y. ))} 
	
	if (discretizeRates) { globalParms <-  .start.discrete.rates(fgyResolution, maxSampleTime = max(bdt$sampleTimes), globalParms ) }
	
	# bdt : binaryDatedTree instance
	if (ncol( bdt$sampleStates ) ==1) 
	  tree <- .calculate.internal.states.unstructuredModel(bdt, maxHeight=maxHeight, globalParms = globalParms)
	else
	  tree_hAs <- .calculate.internal.states(bdt, maxHeight=maxHeight, globalParms = globalParms)
	  tree <- tree_hAs[[1]]
	i<- (length(tree$sampleTimes)+1):(tree$Nnode + length(tree$tip.label))
	if (maxHeight) { 
		internalHeights <- tree$heights[(length(tree$tip.label)+1):length(tree$heights)]
		i <- i[internalHeights <= maxHeight] }
	ll <- sum( log(tree$coalescentRates[i]) ) + sum( tree$logCoalescentSurvivalProbability[i] )#sum( log(tree$coalescentSurvivalProbability[i]) ) 
	if (maxHeight) { ll<- tree$Nnode *  ll/length(i)}
	if (is.nan(ll) | is.na(ll) ) ll <- -Inf
	if (discretizeRates) { globalParms <- .end.discrete.rates(globalParms)}
	#print(paste(date(), 'finish likelihood', ll))
	return (ll)
#~ 	return(list(ll, tree_hAs[[2]], tree$heights[i], log(tree$coalescentRates[i]), log(tree$coalescentSurvivalProbability[i]), tree$logCoalescentSurvivalProbability[i] )) #debug
}

#' Simulate binary dated tree
#' @export
simulatedBinaryDatedTree <- function( x, ...) UseMethod("simulatedBinaryDatedTree")
simulatedBinaryDatedTree.default <- function(sampleTime, sampleStates, FGY=NULL, discretizeRates=FALSE, fgyResolution = 100) 
{
	require(ape)
#~ simulates a coalescent tree, assumes F., G. and Y. are defined
#~ same attributes are defined as binaryDatedTree
	globalParms <- list( USE_DISCRETE_FGY = discretizeRates  )
	if (!is.null(FGY)) { globalParms <- modifyList( globalParms, FGY) }
	else {globalParms <- modifyList( globalParms, list( F. = F., G. = G., Y. = Y. ))} 
#~ <preliminaries>
	n 			<- nrow(sampleStates)
	sampleTimes <- rep( sampleTime, n) 
	Nnode 		<-  n-1
	
	FGY_INDEX <- 1
	
	# NOTE when discretizing rates, this will neglect any changes in rates for t < 0
	if (discretizeRates) 
	{ 
		globalParms <- .start.discrete.rates(fgyResolution, maxSampleTime = max(sampleTimes), globalParms = globalParms) 
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
	calculate.rates <- function(h, parms, globalParms)
	{
		# eliminate diag elements for migration
		t 					<- parms$maxSampleTime - h
		if (discretizeRates)
		{
			.G 					<- globalParms$G.(FGY_INDEX) 
			.F 					<- globalParms$F.(FGY_INDEX)
			.Y 					<- globalParms$Y.(FGY_INDEX) 
		} else
		{
			.G 					<- globalParms$G.(t) 
			.F 					<- globalParms$F.(t)
			.Y 					<- globalParms$Y.(t) 
		}
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
		return( list( X1=X1, X2=X2, lambdaCoalescent=lambdaCoalescent, lambdaMigration=lambdaMigration, lambdaInvisibleTransmission=lambdaInvisibleTransmission) )
	}
	dtheta <- function(h, theta, parms, globalParms, ...){
		if (is.nan(theta)) return(list(NaN))
		if (theta <= parms$r) return(list(NaN)) #terminate early if survival time found
		rates <- calculate.rates(h, parms, globalParms)
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
		rates <- calculate.rates( h, parms, globalParms) 
		parms$rates <- rates
		lambda <- (sum(rates$lambdaCoalescent) + sum(rates$lambdaMigration) + sum(rates$lambdaInvisibleTransmission))
		
		# solve survival time
		# NOTE more accurate to interpolate approxfun( thetaTimes, o[,2] ), & invert to find o[o[,2]==r,1]
		if (discretizeRates){
			eventTime <- rexp(1, rate=lambda) + h
			while (FGY_INDEX < globalParms$FGY_RESOLUTION &  eventTime > globalParms$FGY_H_BOUNDARIES[FGY_INDEX]) {
				if (FGY_INDEX < globalParms$FGY_RESOLUTION &  eventTime > globalParms$FGY_H_BOUNDARIES[FGY_INDEX])  {
					eventTime <- globalParms$FGY_H_BOUNDARIES[FGY_INDEX] 
					FGY_INDEX <- FGY_INDEX+1 ; 
					tryCatch({
								rates <- calculate.rates( eventTime, parms, globalParms) ; parms$rates <- rates
							}, error = function(e) browser() )
					lambda <- (sum(rates$lambdaCoalescent) + sum(rates$lambdaMigration) + sum(rates$lambdaInvisibleTransmission))
					eventTime <- eventTime + rexp(1, rate=lambda)
				}
			}
		} else{
			hub <- h - log(1e-6)/lambda
			thetaTimes <- seq(h, hub, length.out = SIMULATIONTIMERESOLUTION)
			tryCatch( { o <- ode(y=c(theta0), times=thetaTimes, func=dtheta, parms = parms, globalParms = globalParms)}, error =  function(e) {browser() } )
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
		#cumulativeRatesVector
		#indexEvent
		
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
	if (discretizeRates) { globalParms <-  .end.discrete.rates(globalParms)}
	#~ assemble class
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
calculate.cluster.size.moments.from.model <- function(sampleTime, sampleStates, FGY = NULL , maxTime=NA, minTime = NA, timeResolution = 50, discretizeRates=FALSE, fgyResolution = 100 , integrationMethod = 'adams')
{
	require(deSolve)
	globalParms <- list( USE_DISCRETE_FGY = discretizeRates ,  INTEGRATIONMETHOD=integrationMethod )
	if (!is.null(FGY)) { globalParms <- modifyList( globalParms, FGY) }
	else {globalParms <- modifyList( globalParms, list( F. = F., G. = G., Y. = Y. ))} 
	
	n 		<- nrow(sampleStates) 
	m 		<- ncol(sampleStates)
	sampleTimes <- rep(sampleTime, n)
	if (is.na(maxTime)) maxTime <- sampleTime 
	if (is.na(minTime)) minTime <- 0
	heights <- seq( 0, maxTime-minTime, length.out = timeResolution)
	

	if (discretizeRates) 
	{ 
		globalParms <- .start.discrete.rates(fgyResolution, maxSampleTime = max(sampleTimes), globalParms = globalParms) 
	}
	
	
	FGY_INDEX <- 1 
	
	get.fgy <- function(h,t=NULL)
	{
		if (globalParms$USE_DISCRETE_FGY)
		{
			FGY_INDEX <- min( 1+floor( globalParms$FGY_RESOLUTION * h/ globalParms$maxHeight ), globalParms$FGY_RESOLUTION)
			.Y		<- globalParms$Y.(FGY_INDEX)
			.F		<- globalParms$F.(FGY_INDEX)
			.G		<- globalParms$G.(FGY_INDEX)
		} else
		{
			.Y		<- globalParms$Y.(t)
			.F		<- globalParms$F.(t)
			.G		<- globalParms$G.(t)
		}
		list(.F=.F, .G=.G, .Y=.Y,  FGY_INDEX=FGY_INDEX)
	}
	
	dXiXjXij <- function(h, XiXjXij, parms, ...)
	{
		#X are vectors indexed by ancestor type
		Xi	<- XiXjXij[1:m]
		Xj 	<- XiXjXij[(m+1):(2*m)]
		Xij	<- XiXjXij[(2*m+1):(length(XiXjXij))]
		t 	<- parms$maxTime - h
		with(get.fgy(h, t), 
		{
			Xi_Y 	<- Xi / .Y
			Xj_Y 	<- Xj / .Y
			Xij_Y 	<- Xij / .Y
			csFpG 	<- colSums( .F + .G )
			dXi 	<- .F %*% Xi_Y + .G %*% Xi_Y - csFpG * Xi_Y
			dXj 	<- .F %*% Xj_Y + .G %*% Xj_Y - csFpG * Xj_Y
			dXij 	<- .F %*% Xij_Y + .G %*% Xij_Y - csFpG * Xij_Y  + (.F %*% Xj_Y )*Xi_Y  + (.F %*% Xi_Y ) * Xj_Y
			list( c( dXi, dXj, dXij ) )
		})
	}
	
	dA <- function(h, A, parms, ...)
	{
		with(get.fgy(h), 
		{ 
			A_Y 	<- A / .Y
			csFpG 	<- colSums( .F + .G )
			list( .G %*% A_Y - csFpG * A_Y + (.F %*% A_Y) * pmax(1-A_Y, 0) )
		})
	}
	
	
	# solve for A
	A0						<- colSums(sampleStates) 
	parameters 				<- list(maxTime = maxTime)
	o 						<- ode(y=A0, times=heights, func=dA, parms=parameters,  method=globalParms$INTEGRATIONMETHOD) 
	FGY_INDEX 	<- 1
	A 						<- o[,2:(ncol(o))]
	
	# solve Xi Xj and Xij at same time, 3m variables, avoids approxfun nastiness
	# solve second aggregated moment X2 & moments
	Mij_h 		<- array(0, dim = c(m, m, timeResolution), dimnames=list(paste('state',1:m,sep='.'),paste('state',1:m,sep='.'),c()) )
	Mi_h		<- matrix(0, m, timeResolution,dimnames=list(paste('state',1:m,sep='.'),c()))
	for (i in 1:m)
	{ # first tip type
		for (j in i:m)
		{ #second tip type
			#if (i!=j){
			parameters 		<- list(maxTime = maxTime)#, globalParms = globalParms) 
			Xij_h0 			<- rep(0, m)
			if(i==j) 
				Xij_h0[i]	<- A0[i]
			XiXjXij_h0 		<- c( diag(m)[i,] * A0, diag(m)[j,] * A0, Xij_h0) # each X is vector m(ancestor type) X 1
			o 				<- ode(y=XiXjXij_h0, times=heights, func=dXiXjXij,parms=parameters,  method=globalParms$INTEGRATIONMETHOD)
			# moments
			FGY_INDEX 	<- 1
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
	
	if (discretizeRates) { globalParms <- .end.discrete.rates(globalParms)}
	list( heights=heights, Mi= Mi_h, Mij= Mij_h, A = A)
}

calculate.cluster.size.distr.from.tree<- function(bdt, heights)
{
	require(data.table)
	m 		<- ncol( bdt$sampleStates )	
	is.tip 	<- function(bdt, u) 
	{ 
		ifelse(is.na(bdt$daughters[u,1]), TRUE, FALSE)
	}
	identify.tips <- function(bdt, u) 
	{
		if (is.tip(bdt, u)) return(u)
		uv	<- bdt$daughters[u,]
		c(identify.tips(bdt, uv[1]), identify.tips(bdt, uv[2]) ) 
	}
	calculate.sizes.u <- function(bdt, u) 
	{
		tips	<- identify.tips(bdt, u)
		colSums( bdt$sampleStates[tips,,drop=FALSE] )									
	}
	calculate.sizes <- function(bdt, extant)
	{
		if(length(extant)==0)	return( matrix(NA,ncol( bdt$sampleStates ),0) )
		sapply( extant, function(u) calculate.sizes.u(bdt, u) ) 
	}		
	distr	<- lapply( heights, function(h)
			{
				extant 				<- .extant.at.height(h, bdt) 
				sizes_h				<- calculate.sizes(bdt, extant)
				rownames(sizes_h)	<- paste('state',1:m,sep='.')				
				as.data.table( t( rbind(sizes_h, height=ifelse(ncol(sizes_h),h,numeric(0)), ntip=ifelse(ncol(sizes_h),colSums(sizes_h),numeric(0))) ) )
			})
	#print(sapply(distr, function(x) ncol(x) ))
	#print(tail(distr,1))
	distr	<- do.call("rbind",distr)
	distr	
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
	calculate.moment1 <- function(u, i) 
	{
		tips	<- identify.tips(u)
		sum( bdt$sampleStates[tips,i] )				
	}	
	calculate.mean.moment1 <- function(extant, i)
	{
		if(length(extant)==0)	return(NA)
		mean( sapply( extant, function(u) calculate.moment1(u, i) ) ) 
	}
	calculate.moment2 <- function(u, i, j) 
	{
		tips	<- identify.tips(u)
		itips 	<- sum( bdt$sampleStates[tips,i] )
		jtips 	<- sum( bdt$sampleStates[tips,j] )
		itips * jtips
	}
	calculate.mean.moment2 <- function(extant, i, j)
	{
		if (length(extant)==0) return(NA)
		mean( sapply( extant, function(u) calculate.moment2(u, i, j) ) ) 
	}
	for(ih in 1:length(heights))
	{
		h 		<- heights[ih]
		extant 	<- .extant.at.height(h, bdt) 
		if(length(extant) <=1)
		{
			for (i in 1:m)
			{
				Mi[i,ih]	<- Mi[i,ih-1]
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
				Mi[i,ih]	<- calculate.mean.moment1(extant,i)
				for (j in i:m)
				{
					Mij[i,j,ih] <- calculate.mean.moment2(extant,i,j)
					Mij[j,i,ih] <- Mij[i,j,ih]
				}
			}
		}
	}
	
	list( Mi=Mi, Mij=Mij )
}
