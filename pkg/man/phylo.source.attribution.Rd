phylo.source.attribution <- function( trees, sampleTimes, f.t, Y.t, maxTMRCA = NULL, res = 1e3, treeErrorTol = 1e-3)

\name{phylo.source.attribution}
\alias{phylo.source.attribution}
\title{
	Compute infector probabilities between all pairs of sampled lineages. 
}
\description{
	For each pair of sampled lineages (u,v), the probability that u is the donor and v is the recipient is computed. This is conditional on 1) a genealogy representing the evolutionary history of the pathogen with branches in calendar time 2) the times that each lineage are sampled and 3) a time series of incidence and the number infected in the epidemic.  
}

\usage{
phylo.source.attribution( trees, sampleTimes, f.t, Y.t, maxTMRCA = NULL, res = 1e3, treeErrorTol = 1e-3)
}

\arguments{
	\item{trees}{A list of ape::phylo or a multiPhylo}
	\item{sampleTimes}{A named numeric vector providing the time of sampling for each lineage. Should have the same names as trees[[1]]$	tip.label}
	\item{f.t}{ Either a constant providing constant incidence rate of epidemic or a func(t) giving scalar incidence at time t. }
	\item{Y.t}{ Either a constant providing constant number of infected hosts rate of epidemic or a func(t) giving scalar number infected at time t. }
	\item{maxTMRCA}{ Infector probabilities are only computed for pairs with a TMRCA less than this threshold. }
	\item{res}{ Number of time points used for discrete time approximation }
	\item{treeErrorTol}{ Margin of error when computing heights of each node in tree given sampleTimes }
}
\details{
   See vignettes for  detailed usage. 
}
\value{
List with components $W: n by n matrix describing the probability that row infected column
$missingInfector: Length n vector describing probability that infector of each lineage is not sampled
}
\references{
Volz, Erik M., and Simon DW Frost. Inferring the source of transmission with phylogenetic data. (2013): e1003397.
}
\author{
Erik M Volz
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
% \keyword{  }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
