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
mwt(object, grp, log.it = FALSE, localfdr = TRUE)
}
\arguments{
  \item{object}{ Expression data. Either an object of class 'matrix' or 'ExpressionSet'}
  \item{grp}{Numeric or factor. Group indicator variable }
  \item{log.it}{Logical. Take log2 of the data prior to analysis}
  \item{localfdr}{Logical. Shall the function return local FDR (fdr)? Defaults to TRUE}
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
  \item{fdr}{Estimated local FDR from the pvalues}
  \item{df}{degrees of freedom of the MWT test (using central t dist)}
  \item{se2.m}{Moderated standard error}
  \item{d0.prior}{Estimated d0 parameter}
  \item{s2.prior}{Estimated scale parameter for the standard errors}
  \item{lev.stat}{Levene's test statistic}
  \item{lev.FDR}{Levene's test FDR}
}
\references{

  Demissie M, Mascialino B, Calza S, Pawitan Y. Unequal group variances in
  microarray data analyses. Bioinformatics. 2008 May 1;24(9):1168-74.
  PMID: 18344518.

  Ploner A, Calza S, Gusnanto A, Pawitan Y. Multidimensional local false
  discovery rate for microarray studies. Bioinformatics. 2006 Mar
  1;22(5):556-65. PMID: 16368770.
}
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

### Local FDR ###

## Simulate data based on G.Smyth model
require(OCplus)
xdat = MAsim.smyth(ng=10000, p0=0.8, n1=3,n2=3)  ## using Smyth model
dim(xdat)
grp <- factor(colnames(xdat))
colnames(xdat) <- paste("S",1:ncol(xdat),sep=".")

# straight run
out = mwt(xdat, grp)  # get global FDR and local fdr from MWT

## local fdr behaves like fdr2d:
stat = tstatistics(xdat, grp, logse=TRUE)
plot(stat$tstat, stat$logse)
pick = out$fdr<0.1
points(stat$tstat[pick], stat$logse[pick], col='red', pch=16)


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{}
\keyword{}% __ONLY ONE__ keyword per line
