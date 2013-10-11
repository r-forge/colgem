\name{rcolgem-package}
\alias{rcolgem-package}
\alias{rcolgem}
\docType{package}
\title{
Compute the likelihood of a gene genealogy conditonal on a compartmental demographic history.
}
\description{
This package provides methods to compute the likelihood of a gene genealogy using structured coalescent models for populations with complex demographic histories. An example of a complex demographic history is the Susceptible-Infected-Recovered model of infectious disease dynamics and variations thereof. 
}
\details{
\tabular{ll}{
Package: \tab rcolgem\cr
Type: \tab Package\cr
Version: \tab 0.0.1\cr
Date: \tab 2013-10-09\cr
License: \tab GPL-3 \cr
}

\emph{coalescent.log.likelihood}  calculates the coalescent likelihood. One may then proceed to fit the model by maximum likelihood or Bayesian techniques (see, for example, stats4:mle). Please refer to the vignettes for examples of puting a demographic model in the canonical form.

\emph{coalescent.log.likelihood} requires several inputs: 
\itemize{
\item A binary tree with dated tips. \emph{rcolgem} defines the datedBinaryTree class which extends ape:phylo. Instantiating this class requires
	\itemize{
	\item An instance of ape:phylo
	\item  A vector specifying the time of sampling (calendar time) for each taxon in the same order as phylo$tip.label
	\item A matrix specifying the state of each taxon in the same order as phylo$tip.label
	}
\item A compartmental model \eqn{ \mathcal{M} = (F(t), G(t), Y(t))} which is decomposed into the following processes
	\itemize{
	\item $F(t)$ is an \eqn{m\times m} matrix-valued function of time which specifies the population-level rate that gene copies of type $i$ (row) beget gene copies of type $j$ (column)
	\item $G(t)$ is an \eqn{m\times m} matrix-valued function of time which specifies the population-level rate that gene copies of type $i$ (row) become or `migrate to' type $j$ (column)
	\item $Y(t)$ is an m-vector valued function of time that specifies the size of each population
	}
}
The model is specified by assigning the following global variables: 
\code{
F. <<- function(t, ...); 
G. <<- function(t, ...); 
Y. <<- function(t, ...);
}
Once the global variables are defined, and once the dated tree is generated, the likelihood is calculated using 
\code{ coalescent.log.likelihood(tree)  }

}

\author{
Erik M Volz

Maintainer: Erik M Volz <erik.volz+rcolgem@gmail.com
}
\references{
E. M. Volz, Complex population dynamics and the coalescent under neutrality, Genetics, January, 2012
}

\keyword{ package }

\examples{
# load model functions, variables, and MLE fit for a simple HIV model
# Refer to the HIVModel vignette for more details
data(HIVModel)

require(deSolve)

# view model function definitions
F.
G. 
Y. 

# calculate the likelihood of the tree
coalescent.log.likelihood(bdt) 
}
