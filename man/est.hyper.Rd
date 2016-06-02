\name{est.hyper}
\alias{est.hyper}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
est.hyper(z, D, d12)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{z}{ ~~Describe \code{z} here~~ }
  \item{D}{ ~~Describe \code{D} here~~ }
  \item{d12}{ ~~Describe \code{d12} here~~ }
}
\details{
  ~~ If necessary, more details than the description above ~~
}
\value{
  ~Describe the value returned
  If it is a LIST, use
  \item{comp1 }{Description of 'comp1'}
  \item{comp2 }{Description of 'comp2'}
  ...
}
\references{ ~put references to the literature/web site here ~ }
\author{ ~~who you are~~ }
\note{ ~~further notes~~ 

 ~Make other sections like Warning with \section{Warning }{....} ~
}
\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function(z,D,d12){
  f2 = function(d0,D) var(z) - trigamma(D/2) - trigamma(d0/2)
  lim = f2(100,D)
  if (lim<0) D= d12
  d0.est = uniroot(f2,c(1,100), D=D)$root
  s2.est = exp(mean(z) -digamma(D/2) +digamma(d0.est/2) - log(d0.est/D) )
  return(list(d0=d0.est, s2=s2.est))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
