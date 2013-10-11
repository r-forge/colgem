\name{binaryDatedTree}
\alias{binaryDatedTree}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Tree class with dated tips and states for each taxon.
}
\description{
This class extends ape:phylo by adding taxon dates and model-dependent states.
Several helper variables are derived, such as the height of each internal node and mappings from each node to ancestral and daughter nodes.
}
\usage{
binaryDatedTree(x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tree}{
    An object of type ape:phylo. This tree must be binary.
  }
  \item{\dots}{
    sampleTimes : vector of same length as tree$tip.label that defines the time of sampling of each taxon. 
    sampleStates : matrix with number of rows equal to length of tree$tip.label. Each row describes the state of the taxon at the time of sampling. The state is model dependent and should correspond to the F., G. and Y. functions. 
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
Returns a binaryDatedTree object.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Erik M. Volz
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ phylogenetics }
