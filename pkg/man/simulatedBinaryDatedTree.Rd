\name{simulatedBinaryDatedTree}
\alias{simulatedBinaryDatedTree}
\title{Simulate a binary dated tree using a coalescent model given a demographic history}
\usage{
  simulatedBinaryDatedTree(sampleTime, sampleStates, FGY=NULL
    , discretizeRates=FALSE, fgyResolution = 100) 
}
\arguments{
	\item{sampleTime}{ Time when sample takes place.}
	\item{sampleStates}{Numeric matrix of sample states for each taxon. Each row is a probability vector, sums to one, describes probability that each taxon is in each deme.}
	\item{FGY}{List with elements F., G., Y., each is a function of time. F.(t) and G.(t) return a matrix that respectively gives the total rate of births and migrations between demes at time t. Y.(t) returns a vector of the population size in each deme at time t. If NULL, will assume F., G., and Y. are global variables.}
	\item{discretizeRates}{If TRUE, a fast discrete approximation will be used for rates of births and migrations. }
	\item{fgyResolution}{If discretizeRates is TRUE, this variable describes number of time steps in discrete approximation. Larger values give better approximation at greater computational cost. }
}
\value{
	A binaryDatedTree, extends ape::phylo and includes heights of all nodes, time before most recent sample when node occurs. 
}
\description{
}

