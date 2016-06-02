\name{pval2FDR}
\alias{pval2FDR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
pval2FDR(pval, lim = 0.7)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pval}{ ~~Describe \code{pval} here~~ }
  \item{lim}{ ~~Describe \code{lim} here~~ }
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

 ~Make other sections like Warning with Warning{....} ~
}
\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function(pval, lim = 0.7)
{
  n= length(pval)
  Fp = rank(pval)/length(pval)
  p0 = sum(pval>lim)/((1-lim)*n)
  p0 = min(p0, 1)
  
  FDRp = p0 * pmin(pval/Fp, 1)
  ord = order(pval)
  FDR.o = FDRp[ord]
  b = rev(cummin(rev(FDR.o)))
  FDR = rep(0, n)
  FDR[ord] = b
  return(FDR)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
