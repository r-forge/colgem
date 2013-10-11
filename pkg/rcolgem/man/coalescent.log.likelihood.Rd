\name{coalescent.log.likelihood}
\alias{coalescent.log.likelihood}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Compute the likelihood of a gene genealogy conditonal on a compartmental demographic history.
}
\description{
This function expects the following global functions to be defined which represent a canonical decomposition of the demographic model into genealogical processes:
\code{
F. <<- function(t, ...) births
G. <<- function(t, ...) migrations
Y. <<- function(t, ...) size
}
The function also expects a binary genealogy with dated tips. The taxa of the tree must also have defined states corresponding to the given demographic history (see \emph{binaryDatedTree}). 

}
\usage{
coalescent.log.likelihood(bdt, integrationMethod = "rk4")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{bdt}{
    A binaryDatedTree object which extends ape:phylo and includes taxon dates and states.
}
  \item{integrationMethod}{
     Passed to lsoda. Higher precision may be obtained with methods such as \emph{adams} at some computational cost. 
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
 The log likelihood, type numeric. 
}
\references{
E. M. Volz, Complex population dynamics and the coalescent under neutrality, Genetics, January, 2012
}
\author{
Erik M Volz
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{binaryDatedTree}}
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as...
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ coalescent }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
