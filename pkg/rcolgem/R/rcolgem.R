#October 28 2013
# expects F.(t), G.(t), and Y.(t) to be in namespace
#~ TODO finite size correction for pair (i,j) of lineages at internal node (see written notes): 
#~	if i transmits and is type k: 
#~ 		pjk -> pjk * (Yk-1)/Yk
#~ 		pjl (l\neq k) -> pjl * Yk/(Yk-pjk)
#~ TODO option to correct for direct ancestor sampling if doing serial samples and there is a 0-branch length
#~ 	add these terms to the likelihood:
#~ 		if 0 bl at s_i: \sum_k pik Ak / Yk
#~ 		else: \sum_k pik (1-Ak/Yk)

require(ape)
require(deSolve)

#PRELIMINARIES 
calculate.heights <- function(phylo){
	phylo$maxSampleTime  = phylo$T <- max(phylo$sampleTimes)
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
			nextgen <- c(nextgen, u)
		}
		curgen <- unique(nextgen)
	}
	phylo$heights <- heights
	return(phylo)
}

calculate.edgemap <- function(phylo){
	inEdgeMap <- rep(-1, length((phylo$Nnode + length(phylo$tip.label))))
	outEdgeMap <- matrix(-1, (phylo$Nnode + length(phylo$tip.label)), 2)
	parent <- 1:(phylo$Nnode + length(phylo$tip.label)) 
	daughters <- matrix(-1, (phylo$Nnode + length(phylo$tip.label)), 2)
	for (u in 1:(phylo$Nnode + length(phylo$tip.label))){
		if (u!=length(phylo$tip.label)+1){ #if u not root
			inEdgeMap[u] <- which( phylo$edge[,2]==u ) 
		} else{ 
			inEdgeMap[u] <- u
		}
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

initialize.states <- function(phylo){
	phylo$m =m <- dim(phylo$sampleStates)[2]
	phylo$lstates <- matrix(-1, (phylo$Nnode + length(phylo$tip.label)), m)
	phylo$mstates <- matrix(-1, (phylo$Nnode + length(phylo$tip.label)), m)
	phylo$ustates <- matrix(-1, (phylo$Nnode + length(phylo$tip.label)), m)
	phylo$lstates[1:nrow(phylo$sampleStates),] <- phylo$sampleStates
	phylo$mstates[1:nrow(phylo$sampleStates),] <- phylo$sampleStates
	phylo$ustates[1:nrow(phylo$sampleStates),] <- phylo$sampleStates
	return(phylo)
}

# binaryDatedTree class, includes heights for each node and other helper variables
binaryDatedTree <- function( x, ...) UseMethod("binaryDatedTree") 
binaryDatedTree.default <- function( phylo, sampleTimes, sampleStates){
	if (phylo$Nnode != length(phylo$tip.label) - 1 ) { stop('Object class phylo is not a binary tree.') }
	phylo$sampleTimes <- sampleTimes
	phylo$sampleStates <- sampleStates
	phylo <- calculate.heights(phylo)
	phylo <- calculate.edgemap(phylo)
	phylo <- initialize.states(phylo)
	phylo$coalescentRates <- rep(0, (phylo$Nnode + length(phylo$tip.label)))
	phylo$coalescentSurvivalProbability <- rep(1, (phylo$Nnode + length(phylo$tip.label)))
	class(phylo) <- "binaryDatedTree"
	return(phylo)
}


extant.at.height <- function(h, tree){
{return( which( tree$heights <= h & tree$parentheights > h)  )}
}


##########################
#CALCULATE INTERNAL STATES

dPiAL <- function(h, y, parms, ...){
	# conditional on no coalescent
	p <- y[1:(parms$m-1)]
	A <- y[(parms$m):(2*parms$m-1)]
	L <- y[length(y)] #cumulative hazard
	pp <- c(p, 1 - sum(p))
	t <- parms$treeT - h
	
	.G <- G.(t) 
	.F <- F.(t)
	.Y <- Y.(t) 
	
	X1 <- A / .Y
	X2 <-  (.Y - A ) / .Y
	
	FklXAk_Yk <- (.F * X1)
	if(FINITESIZECORRECTIONS){
		diag(FklXAk_Yk) <- diag(FklXAk_Yk) * .Y / (pmax(.Y,1.01)-1)
		}
	dL <- sum( FklXAk_Yk %*% X1 ) 
	
	#TODO would be faster to solve m- 1 equations since sum(A) is conserved
	dA <- c(- as.vector(t(.G) %*% X1) 
	+ as.vector(.G %*% X1)
	- as.vector( t(.F) %*% X2 ) * X1
	+ as.vector( .F %*% X1) * X2)
	
	R <- t(.G) / .Y + t(.F * X2) / .Y
	R <- R * (1- diag(length(.Y)))
	R <- R  -rowSums(R) * diag(length(.Y))
	
	dPi <-  (t(R) %*% pp)[1:(length(.Y)-1)]
	
	return(list( c(dPi, dA, dL) ))
}


solve.Pi.and.AL <- function(h0, h1, A0, L0, tree){
	P0 <- diag(tree$m)[,1:(tree$m-1)]
	if (tree$m <=2) P0 <- t(t(P0))
	parameters <- list(treeT = tree$T, m = tree$m)
	Pi1s <- c()
	for (i in 1:tree$m){
		y0 <- c( P0[i,], A0, L0)
		out0 <- ode(y=y0, times=c(h0, h1), func=dPiAL, parms=parameters, method=INTEGRATIONMETHOD ) 
		Pi1 <- abs(out0[nrow(out0),2:(1 + tree$m-1)])
		if( sum(Pi1) > 1) {Pi1 <- Pi1 / sum(Pi1)}
		Pi1s <- rbind(Pi1s, Pi1 )
	}
	Pi1s <- cbind(Pi1s, 1 - rowSums(Pi1s))
	A1 <- out0[nrow(out0), (tree$m + 1):(2*tree$m-1 + 1)]
	L1 <- out0[nrow(out0), ncol(out0)]
	return ( list(  unname(Pi1s), unname(A1), unname (L1) ) )
}

calculate.internal.states <- function(tree){
	eventTimes <- unique( sort(tree$heights) )
	S <- 1
	L <- 0
	for (ih in 1:(length(eventTimes)-1)){
		h0 <- eventTimes[ih]
		h1 <- eventTimes[ih+1]
		#get A0, process new samples, calculate lstate for new lines
		extantLines <- extant.at.height(h0, tree)
		if (length(extantLines) > 1 ){
			A0 <- colSums(tree$mstates[extantLines,])
		} else{ A0 <- tree$mstates[extantLines,] }
		
		#solve P[i,], nlft, and S=exp(-L)
		out <- solve.Pi.and.AL(h0, h1, A0, L, tree)
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
		.F <- F.(tree$T - h1)
		.Y <- Y.(tree$T - h1)
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
			ratekl <- FklXpuk_Yk %*% (tree$ustates[v,]/.Y) + FklXpvk_Yk %*% (tree$ustates[u,]/.Y)
			
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
						fterm <- p_i * p_a * (.Y-1)/(.Y-p_a)
						smat <- t( matrix( rep( p_a * .Y/(.Y-p_i), tree$m), nrow=tree$m) )
						diag(smat) <- 0
						sterm <- p_i * rowSums(smat)
						tree$mstates[i,] <- fterm + sterm
						tree$mstates[i,] <- tree$mstates[i,] / sum(tree$mstates[i,])
					}
				}
			}
		}
		
		if (length(newNodes) > 0) {L<-0}
	}
	return(tree)
}


#############################
# CALCULATE LIKELIHOOD
coalescent.log.likelihood <- function(bdt, integrationMethod = 'rk4', finiteSizeCorrections=FALSE){
	# bdt : binaryDatedTree instance
	INTEGRATIONMETHOD <<- integrationMethod
	FINITESIZECORRECTIONS <<- finiteSizeCorrections
	tree <- calculate.internal.states(bdt)
	i<- (length(tree$sampleTimes)+1):(tree$Nnode + length(tree$tip.label))
	ll <- sum( log(tree$coalescentRates[i]) ) + sum( log(tree$coalescentSurvivalProbability[i]) )
	if (is.nan(ll) | is.na(ll) ) ll <- -Inf
	
	return(ll)
}
