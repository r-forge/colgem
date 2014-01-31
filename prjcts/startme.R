#! /Library/Frameworks/R.framework/Versions/2.15/Resources/bin/Rscript
##	the first line must be the absolute path to Rscript - need to set manually
##! /apps/R/2.15/lib64/R/bin/Rscript
###############################################################################
#
#
# Author: oliver ratmann
# file: startme.R
#
# usage from R:
#	> setwd("/Users/Oliver/workspace_sandbox/colgem/prjcts")
#	> source("startme.R")
# usage from bash:
#> startme.R 
#
#	README:	
#	the script DO NOT currently work for some reason from the command line, and only work when started in R.
#	TODO	Erik, the problem might be related to the usage of global functions. 
#			It might be good to fix this.
#	SO:
#	open an R window
#	type setwd("YOURPATH/colgem/prjcts")
#	type source("startme.R")
#	copy and paste the script functions, eg 'cg.sde'
#
###############################################################################
#	set these as needed
CODE.HOME	<<- "/Users/Oliver/workspace_sandbox/colgem"
HOME		<<- "/Users/Oliver/workspace_sandbox/colgem/prjcts"
#CODE.HOME	<<- "/work/or105/libs/colgem"
#HOME		<<- "/work/or105/libs/colgem/prjcts"

default.fun	<- "colgem.project.moments.for.odesystem"
default.fun	<- "cg.sde"
default.fun	<- "cg.pipeline"
###############################################################################
args <- commandArgs()
if(!any(args=='--args'))
	args<- vector("numeric",0)
if(any(args=='--args'))
	args<- args[-(1:match("--args", args)) ]
###############################################################################
DATA		<<- HOME
SCRIPT		<<- HOME
DEBUG		<<- 0
LIB.LOC		<<- NULL
#LIB.LOC		<<- paste(CODE.HOME,"../",sep='')
EPS			<<- 1e-12
###############################################################################
require(rcolgem)
###############################################################################
#	load all R functions and project source files. This avoids to rebuild the package all the time.
function.list	<-c(	list.files(path= paste(CODE.HOME,"pkg/R",sep='/'), pattern = ".R$", all.files = FALSE, full.names = TRUE, recursive = FALSE),
						list.files(path= paste(HOME,"misc",sep='/'), pattern = ".R$", all.files = FALSE, full.names = TRUE, recursive = FALSE)
						)
sapply(function.list,function(x) source(x,echo=FALSE,print.eval=FALSE, verbose=FALSE))
###############################################################################
my.mkdir<-function(root,data.name)
{
	if(length(dir(root,pattern=paste('^',data.name,'$',sep='')))==0)
		system(paste("mkdir ",paste(root,data.name,sep='/'),sep=''))
}
###############################################################################
my.dumpframes<- function()
{
	geterrmessage()
	dump.frames()
	cat(paste("\nmy.dumpframes dump 'last.dump' to file",paste(DATA,paste("debug_",paste(strsplit(date(),' ')[[1]],collapse='_'),".rda\n",sep=''),sep='/')))
	save(last.dump, file=paste(DATA,paste("debug_",paste(strsplit(date(),' ')[[1]],collapse='_'),".rda",sep=''),sep='/'))
	q()
}
###############################################################################
my.make.documentation<- function()
{
	require(roxygen2)		
	roxygenize(CODE.HOME)
}
###############################################################################
my.fade.col<-function(col,alpha=0.5)
{
	return(rgb(col2rgb(col)[1]/255,col2rgb(col)[2]/255,col2rgb(col)[3]/255,alpha))
}
###############################################################################
print.v<- function(x,cut=3,digits=4,prefix= "simu_",print.char= TRUE, as.R= TRUE)
{
	if(as.R)
	{
		tmp<- paste("c(",paste(c(x,recursive=T),collapse=',',sep=''),')',sep='')
		if(!is.null(names(x)))
			tmp<- paste("{tmp<-", tmp, "; names(tmp)<- ", paste('c("',paste(c(names(x),recursive=T),collapse='", "',sep=''),'")',sep=''), "; tmp}", sep= '', collapse= '')
	}
	else
	{
		if(!is.null(names(x)))
		{
			m<- matrix(NA,nrow=2,ncol=length(x))
			m[1,]<- substr(names(x),1,cut)
			m[2,]<- round( x, digits=digits )
			if(cut==0)		m<- m[2,]
			tmp<- gsub('.',',',paste(prefix,paste(as.vector(m), collapse='_',sep=''),sep=''),fixed=T)
		}
		else
			tmp<- gsub('.',',',paste(prefix,paste(round( x, digits=digits ), collapse='_',sep=''),sep=''),fixed=T)
	}
	if(print.char) print(tmp)
	tmp
}
###############################################################################
plot.2D.dens<- function(x,y,xlab,ylab,xlim=NA,ylim=NA,nbin=NA,width.infl=2,n.hists=5,method="gauss", palette= "topo", persp.theta= -30, persp.phi= 30, zero.abline=TRUE, ...)
{
	if(!method%in%c("gauss","ash","persp"))	stop("plot.2D.dens: exception 1a")
	if(!palette%in%c("topo","heat","gray"))	stop("plot.2D.dens: exception 1b")
	switch(method,
			gauss=
					{
						require(KernSmooth)
						require(fields)
						x.bw<- width.infl*diff(summary(x)[c(2,5)])
						y.bw<- width.infl*diff(summary(y)[c(2,5)])
						if(!x.bw) x.bw<- EPS
						if(!y.bw) y.bw<- EPS
						if(any(is.na(xlim)))	xlim<- range(x)+c(-1.5,1.5)*x.bw
						if(any(is.na(ylim)))	ylim<- range(y)+c(-1.5,1.5)*y.bw
						dens <- bkde2D(cbind(x, y), range.x=list(xlim,ylim),bandwidth=c(x.bw,y.bw))
						contour(dens$x1, dens$x2, dens$fhat,xlab=xlab,ylab=ylab)
						if(zero.abline) abline(v=0,col="black",lty=3,lwd=1.5)
						if(zero.abline) abline(h=0,col="black",lty=3,lwd=1.5)
					},
			ash={
				require(ash)
				if(any(is.na(xlim))) xlim<- range(x)*1.05
				if(any(is.na(ylim))) ylim<- range(y)*1.05
				if(any(is.na(nbin))) nbin<- 2*c(nclass.Sturges(x),nclass.Sturges(y))
				bins<- bin2(cbind(x, y), ab=rbind(xlim,ylim),nbin=nbin)
				f <- ash2(bins,rep(n.hists,2))
				#image(f$x,f$y,f$z, col=rainbow(50,start= 4/6,end=0),xlab=xlab,ylab=ylab)
				if(palette=="topo")					image(f$x,f$y,f$z, col=tail( topo.colors(trunc(50*1.4)), 50 ),xlab=xlab,ylab=ylab,...)
				else if(palette=="gray")			image(f$x,f$y,f$z, col=head( rev(gray(seq(0,.95,len=trunc(50*1.4)))), 50),xlab=xlab,ylab=ylab,...)					
				else								image(f$x,f$y,f$z, col=heat.colors( 50 ),xlab=xlab,ylab=ylab,...)
				contour(f$x,f$y,f$z,add=TRUE, nlevels= 5)
				if(zero.abline) abline(v=0,col="black",lty=2,lwd=2)
				if(zero.abline) abline(h=0,col="black",lty=2,lwd=2)
			},
			persp={
				require(ash)
				require(MASS)
				if(any(is.na(xlim))) xlim<- range(x)*1.05
				if(any(is.na(ylim))) ylim<- range(y)*1.05
				bins<- bin2(cbind(x, y), ab=rbind(xlim,ylim),nbin=2*c(nclass.Sturges(x),nclass.Sturges(y)))
				f <- ash2(bins,rep(n.hists,2))
				
				nrz <- nrow(f$z)
				ncz <- ncol(f$z)
				col<- tail(topo.colors(trunc(1.4 * 50)),50)
				fcol      <- col[trunc(f$z / max(f$z)*(50-1))+1]
				dim(fcol) <- c(nrz,ncz)
				fcol      <- fcol[-nrz,-ncz]
				par(mar=c(1/2,1/2,1/2,1/2))
				persp(x=f$x,y=f$y,z=f$z,col= fcol,theta=persp.theta,phi=persp.phi,xlab=xlab,ylab=ylab,zlab='', ticktype= "detailed" )
			})
}
###############################################################################
plot.persplocfit<- function(x, pv, theta= 30, phi= 20, palette= "gray",tcl=-0.05,...)
{	
	d <- x$mi["d"]
	ev <- x$mi["ev"]
	where <- "grid"
	
	pv <- match(pv, x$vnames)
	tv <- (1:d)[-pv]
	vrs <- c(pv, tv)
	if (any(duplicated(vrs))) 
		warning("Duplicated variables in pv, tv")
	if (any((vrs <= 0) | (vrs > d))) 
		stop("Invalid variable numbers in pv, tv")
	m <- ifelse(d == 1, 100, 40) 
	m <- rep(m, d)
	m[tv] <- mtv<- 6		
	xl <- x$box
	marg <- lfmarg(xl, m)
	pred <- locfit:::preplot.locfit(x, marg, band = "none", tr = NULL, what = "coef", get.data = 0, f3d = 0)
	z<- matrix(pred$fit, nrow = length(marg[[1]]))
	nbcol <- 100
	if(palette=="gray")	
		color<- head( rev(gray(seq(0,0.95,len=trunc(nbcol*1.4)))), nbcol)
	else
	{
		jet.colors <- colorRampPalette( c("blue", "green") )
		color <- jet.colors(nbcol)
	}
	# Compute the z-value at the facet centres
	nrz <- nrow(z)
	ncz <- ncol(z)
	zfacet <- z[-1, -1] + z[-1, -ncz] + z[-nrz, -1] + z[-nrz, -ncz]
	# Recode facet z-values into color indices
	facetcol <- cut(zfacet, nbcol)		
	par(mar=c(0,2.5,0,0), tcl=tcl)
	pmat<- persp(marg[[1]], marg[[2]], z, zlim= range(z)*1.1, col = color[facetcol], shade= 0.1, border=NA, ticktype = "detailed", ltheta = 120,theta = theta, phi = phi, expand = 0.75, box=1, ... )
	
	list(pmat=pmat, x= marg[[1]], y= marg[[2]], z= z)
}
###############################################################################
#
#
#
###############################################################################
argv<- list()
if(length(args))
{
	tmp<- na.omit(sapply(args,function(arg)
					{
						switch(substr(arg,2,4),
								exe= return(substr(arg,6,nchar(arg))),
								NA)
					}))
	if(length(tmp)!=0)
	{
		if(length(tmp)>1) stop("abc-n.startme.R: duplicate -exe")
		else default.fun<- switch(tmp[1],
					MAKE.DOCUMENTATION		 = "my.make.documentation",
					LKL.SDE					 = "cg.sde.fulllkl",
					MOM.SDE.PSEUDODATA		 = "cg.sde.get.pseudodata",
					MOM.ODE					 = "colgem.project.moments.for.odesystem",
					MOM.SDE					 = "cg.sde",					
					MOM.SDE.SIMMO			 = "cg.sde.get.mM",
					MOM.SDE.LN3D			 = "cg.sde.eM.pseudolkl.fit3D"					
					)
	}
	tmp<- na.omit(sapply(args,function(arg)
					{
						switch(substr(arg,2,6),
								debug= 1,
								NA)
					}))		
	if(length(tmp)!=0)	DEBUG<<- tmp[1]	
	argv<<- args
}
###############################################################################
if(DEBUG)	options(error= my.dumpframes)	
cat(paste("\ncolgem.startme.R: ",ifelse(DEBUG,"debug",""),"call",default.fun))
do.call(default.fun,list()) 	
cat("\ncolgem: ",ifelse(DEBUG,"debug","")," end\n")
