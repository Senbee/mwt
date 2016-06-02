\name{MAsim.uneqvar}
\alias{MAsim.uneqvar}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
MAsim.uneqvar(ng = 10000, n = 3, n1 = n, n2 = n, p0 = 0.9, d01 = 10, s2_01 = 1, v0m = 1, v0var = 1, p0var = 0.7)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ng}{ ~~Describe \code{ng} here~~ }
  \item{n}{ ~~Describe \code{n} here~~ }
  \item{n1}{ ~~Describe \code{n1} here~~ }
  \item{n2}{ ~~Describe \code{n2} here~~ }
  \item{p0}{ ~~Describe \code{p0} here~~ }
  \item{d01}{ ~~Describe \code{d01} here~~ }
  \item{s2_01}{ ~~Describe \code{s2_01} here~~ }
  \item{v0m}{ ~~Describe \code{v0m} here~~ }
  \item{v0var}{ ~~Describe \code{v0var} here~~ }
  \item{p0var}{ ~~Describe \code{p0var} here~~ }
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
function (ng = 10000, n = 3, n1 = n, n2 = n, p0 = 0.9, 
                           d01 = 10, s2_01 = 1, v0m = 1, v0var=1, p0var=0.7) 
{
  nn = n1 + n2
  group = rep(c(0, 1), c(n1, n2))
  s2 = d01*s2_01/rchisq(ng, df = d01)
  s2_g1 = s2*exp(rnorm(ng,mean=0,sd= sqrt(v0var)))
  s2_g2 = s2*exp(rnorm(ng,mean=0,sd= sqrt(v0var)))
  
  eqvar = runif(ng)< p0var  # equal variance
  s2_g1[eqvar] = s2[eqvar]
  s2_g2[eqvar] = s2[eqvar]
  s2_g = (n1*s2_g1 + n2*s2_g2)/(n1+n2)   # average variance
  
  xdat1 = matrix(rnorm(ng*n1, sd = sqrt(s2_g1)), ncol = n1)
  
  ndx = runif(ng) > p0
  nde = sum(ndx)
  xmean2 = rep(0, ng)
  xmean2[ndx] = rnorm(nde, mean = 0, sd = sqrt(v0m * s2_g))
  xdat2 = matrix(rnorm(ng*n2, sd = sqrt(s2_g2)), ncol = n2) + xmean2
  
  xdat = cbind(xdat1, xdat2)
  colnames(xdat) = as.character(group)
  
  des = des.var = rep(FALSE, ng)   # DE status for mean and variance
  des[ndx] = TRUE
  des.var[!eqvar] = TRUE
  attr(xdat, "DE") = des
  attr(xdat, "DE.var") = des.var
  
  xdat
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
