\name{mwt}
\alias{mwt}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Function to compute global FDR based on the moderated Welch test }
\description{
  MWT compares two independent groups using t-test. It is designed
  to deal with small-sample microarray data where the group
  variances might be unequal. In small samples it is better than
  either the standard t-test or its Welch version.
}
\usage{
mwt(xdat, grp, log.it = FALSE)
}
\arguments{
  \item{xdat}{ Expression data. Either an object of class 'matrix' or 'ExpressionSet'}
  \item{grp}{Numeric or factor. Group indicator variable }
  \item{log.it}{Logical. Take log2 of the data prior to analysis}
 }
\details{
  The statistic is equal mean difference divided by its standard error,
  where the std error is a weighted average of pooled and unpooled 
  standard errors, and the weight is the FDR for equal variance. The 
  std error is further penalized to avoid small values.
}
\value{
A list containing
  \item{MWT}{Moderated Welch statistic}
  \item{pvalue}{Corresponding p-values from MWT}
  \item{FDR}{Estimated global FDR from the pvalues}
  \item{df}{degrees of freedom of the MWT test (using central t dist)}
  \item{s2.wm}{Posterior estimate of variance}
  \item{d0.prior}{Estimated d0 parameter}
  \item{s2.prior}{Estimated scale parameter for the standard errors}
  \item{lev.stat}{Levene's test statistic}
  \item{lev.FDR}{Levene's test FDR}
}
\references{
  Demissie M, et al: Effects of unequal variance in microarray data analysis }
\author{Pawitan Y and Calza S}
\note{ }
\seealso{ }
\examples{
# simulate data with unequal variance
xdat = MAsim.uneqvar(ng=10000,n1=3,n2=9)  
dim(xdat)
grp <- factor(colnames(xdat))
colnames(xdat) <- paste("S",1:ncol(xdat),sep=".")

# straight run
out = mwt(xdat, grp)  # get FDR from MWT
names(out)
plot(out$MWT, out$FDR)

# alternative run using ExpressionSet class
eset <- new("ExpressionSet",exprs=xdat,
            phenoData=new("AnnotatedDataFrame",
              data=data.frame(GRP=grp,row.names=colnames(xdat))))
out = mwt(eset, "GRP")  # get FDR from MWT
plot(out$MWT, out$FDR)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{}
\keyword{}% __ONLY ONE__ keyword per line
