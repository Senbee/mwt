
## -- finding d0 and s0

## est.hyper <- function(z,D,d12)
## {
##   f2 <- function(d0,D)
##     {
##       var(z) - trigamma(D/2) - trigamma(d0/2)
##     }
  
##   lim <- f2(100,D)

##   if(lim < 0)
##     D <- d12

##   d0.est = uniroot(f2,c(1,100), D=D)$root
##   s2.est = exp(mean(z) -digamma(D/2) +digamma(d0.est/2) - log(d0.est/D) )
  
##   return(list(d0=d0.est, s2=s2.est))
## }

## allows for NA
est.hyper <- function (z, D, d12) 
{
  f2 <- function(d0, D)
    {
      var(z,na.rm=TRUE) - trigamma(D/2) - trigamma(d0/2)
    }

  lim = f2(100, D)

  if (lim<0)
    d0.est <- 100   # f2(1) and f2(100) have same sign

  if (lim>0)
    d0.est <- uniroot(f2, c(1, 100), D = D)$root

  s2.est <- exp(mean(z,na.rm=TRUE) - digamma(D/2) + digamma(d0.est/2) - log(d0.est/D))

  return(list(d0 = d0.est, s2 = s2.est))
}



## -- transform pval to FDR

## pval2FDR <- function(pval, lim = 0.7)
## {
##   n= length(pval)
##   Fp = rank(pval)/length(pval)
##   p0 = sum(pval>lim)/((1-lim)*n)
##   p0 = min(p0, 1)
  
##   FDRp = p0 * pmin(pval/Fp, 1)
##   ord = order(pval)
##   FDR.o = FDRp[ord]
##   b = rev(cummin(rev(FDR.o)))
##   FDR = rep(0, n)
##   FDR[ord] = b
##   attr(FDR,"p0") <- p0
##   return(FDR)
## }

## allows for NA
pval2FDR <- function (pval, lim = 0.7) 
{
  
    
  n1 = length(pval)
  ok.id <- 1:n1
  
  if(any(is.na(pval)))
    {
      ok.id <- which(!is.na(pval))
      pval <- na.omit(pval)
    }

  n = length(pval)
  Fp = rank(pval)/length(pval)
  p0 = sum(pval > lim)/((1 - lim) * n)
  p0 = min(p0, 1)
  FDRp = p0 * pmin(pval/Fp, 1)
  ord = order(pval)
  FDR.o = FDRp[ord]
  b = rev(cummin(rev(FDR.o)))
  FDR = rep(0, n)
  FDR[ord] = b

  out.FDR <- rep(NA, n1)
  out.FDR[ok.id] <- FDR
  
  return(out.FDR)
}





## simulation of microarray data with unequal variance
## as used in Demissie et al 2007


MAsim.uneqvar <- function (ng = 10000, n = 3, n1 = n, n2 = n, p0 = 0.9, 
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

#
# simulating 3 groups
#
# there is a common variance component across the groups
#
MAsim.var6 <- function (ng = 10000, n = 10, n1 = n, n2 = n, n3=n, p0 = 0.9, 
                        d01 = 4, s2_01 = 4, v0m = 2, v0var=1, p0var=0.9) 
{
  nn = n1 + n2 + n3
  group = rep(c(0, 1, 2), c(n1, n2, n3))
  s2 = d01*s2_01/rchisq(ng, df = d01)
  s2_g1 = s2*exp(rnorm(ng,mean=0,sd= sqrt(v0var)))
  s2_g2 = s2*exp(rnorm(ng,mean=0,sd= sqrt(v0var)))
  s2_g3 = s2*exp(rnorm(ng,mean=0,sd= sqrt(v0var)))
  
  eqvar = runif(ng)< p0var  # equal variance
  s2_g1[eqvar] = s2[eqvar]
  eqvar = runif(ng)< p0var  # equal variance
  s2_g2[eqvar] = s2[eqvar]
  eqvar = runif(ng)< p0var  # equal variance
  s2_g3[eqvar] = s2[eqvar]
  
  s2_g = (n1*s2_g1 + n2*s2_g2 + n3*s2_g3)/(n1+n2+n3)   # average variance
  
  xdat1 = matrix(rnorm(ng*n1, sd = sqrt(s2_g1)), ncol = n1)
  
  ndx2 = runif(ng) > p0
  nde = sum(ndx2)
  xmean2 = rep(0, ng)
  xmean2[ndx2] = rnorm(nde, mean = 0, sd = sqrt(v0m * s2_g))
  xdat2 = matrix(rnorm(ng*n2, sd = sqrt(s2_g2)), ncol = n2) + xmean2
  
  ndx3 = runif(ng) > p0
  nde = sum(ndx3)
  xmean3 = rep(0, ng)
  xmean3[ndx3] = rnorm(nde, mean = 0, sd = sqrt(v0m * s2_g))
  xdat3 = matrix(rnorm(ng*n3, sd = sqrt(s2_g3)), ncol = n3) + xmean3
  
  xdat = cbind(xdat1, xdat2, xdat3)
  colnames(xdat) = as.character(group)
  
  des = des.var = rep(FALSE, ng)   # DE status for mean and variance
  des[(ndx2 | ndx3)] = TRUE
  des.var[!eqvar] = TRUE
  attr(xdat, "DE") = des
  attr(xdat, "DE.var") = des.var
  
  xdat
}

## .............................  true FDP for F test from ANOVA
## stat = F statistic; 
## DE = true DE status

FDP.F <- function(stat,DE)
{   
  ng = length(stat)
  ord = order(-stat)   # to get at the large statistics first
  TDP = cumsum(DE[ord])/c(1:ng)
  FDP = 1- cummin(TDP)
  return(FDP)
}



## --- Levene's test
## levene <- function(xdat,grp)
## {
##   glab = unique(grp)
##   ngr = length(glab)
##   n = mn = s2 = NULL
  
##   ## center each group first
##   X0 = NULL
##   for (i in 1:ngr)
##     { 
##       ndx = grp == glab[i]
##       mni = rowMeans(xdat[, ndx])
##       x0 = xdat[, ndx] - mni
##       X0 = cbind(X0,x0)
##     } 
  
##   ## redefine xdat as the abs deviations
##   xdat = abs(X0)
##   for (i in 1:ngr){ 
##     ndx = grp == glab[i]
##     ni = sum(ndx)
##     mni = rowMeans(xdat[, ndx])
##     x0 = xdat[, ndx] - mni
##     s2i = rowSums(x0 * x0)/(ni-1)
##     n = c(n, ni)
##     mn = cbind(mn, mni)
##     s2 = cbind(s2, s2i)
##   } 
  
##   ## .............   F statistic
##   N = sum(n)
##   mmn = rowSums(xdat)/N
##   mn0 = mn - mmn
##   num = rowSums(t(t(mn0*mn0)*n))/(ngr-1)
##   den = rowSums(t(t(s2*(n-1))))/(N-ngr)
##   F3 = num/den
##   pval = pf(F3, df1= ngr-1, df2 = N-ngr)
##   pval = ifelse(pval<0.5,2*pval, 2*(1-pval))
##   lvn.FDR = pval2FDR(pval)
  
##   return(list(pvalue= pval, FDR= lvn.FDR))
## }

## allows for NA

levene <- function(xdat,grp,na.rm=TRUE)
{
  glab = unique(grp)
  ngr = length(glab)
  n = mn = s2 = NULL
  
  ## center each group first
  X0 = NULL
  for (i in 1:ngr)
    { 
      ndx = grp == glab[i]
      mni = rowMeans(xdat[, ndx],na.rm=na.rm)
      x0 = xdat[, ndx] - mni
      X0 = cbind(X0,x0)
    } 
  
  ## redefine xdat as the abs deviations
  xdat = abs(X0)
  for (i in 1:ngr){ 
    ndx = grp == glab[i]
    ni = sum(ndx)
    mni = rowMeans(xdat[, ndx],na.rm=na.rm)
    x0 = xdat[, ndx] - mni
    s2i = rowSums(x0 * x0, na.rm=na.rm)/(ni-1)
    n = c(n, ni)
    mn = cbind(mn, mni)
    s2 = cbind(s2, s2i)
  } 
  
  ## .............   F statistic
  N = sum(n)
  mmn = rowSums(xdat, na.rm=na.rm)/N
  mn0 = mn - mmn
  num = rowSums(t(t(mn0*mn0)*n))/(ngr-1)
  den = rowSums(t(t(s2*(n-1))))/(N-ngr)
  F3 = num/den
  pval = pf(F3, df1= ngr-1, df2 = N-ngr)
  pval = ifelse(pval<0.5,2*pval, 2*(1-pval))
  lvn.FDR = pval2FDR(pval)
  
  return(list(statistic=F3,pvalue= pval, FDR= lvn.FDR))
}
