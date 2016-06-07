setGeneric("mwt",
           function(object, grp, log.it=FALSE, locfdr = FALSE)
           {
             standardGeneric("mwt")
           })


setGeneric("mwt2",
           function(object, grp, log.it=FALSE)
           {
             standardGeneric("mwt2")
           })


setMethod("mwt",
    signature(object="ExpressionSet"),
    function (object, grp, log.it=FALSE, locfdr = FALSE)
          {
            if(length(grp) == 1)
              if(!is.character(grp))
                stop("grp must be either a character vector or a vector of length = ",
                     ncol(exprs(object)))
              else
                if(!grp %in% names(pData(object)))
                  stop("grp must be an existing variable in object's phenoData")
                else
                  grp <- object[[grp]]

            if(length(unique(grp)) > 2)
              stop("Not a two-groups comparison. Use mwt2 for multiple groups")

            if(log.it)
              xdat <- log2(exprs(object))
            else
              xdat <- exprs(object)
            
            ans <- stat2(xdat,grp,locfdr=locfdr)
            
            return(ans)
          }
          
          )


setMethod("mwt",signature(object="matrix"),
          function(object, grp, log.it=FALSE,locfdr=FALSE)
          {
            if(log.it)
              object <- log2(object)
            
            if(length(grp) != ncol(object))
              stop("grp must be a vector of length ",ncol(object))
            
            if(length(unique(grp)) > 2)
              stop("Not a two-groups comparison. Use mwt2 for multiple groups")
            ans <- stat2(object,grp,locfdr=locfdr)
            
            return(ans)
            
            
          }
          )



setMethod("mwt2",
    signature(object="matrix"),
    function (object, grp, log.it=FALSE) 
          {

                        
            if(length(grp) == 1)
              if(!is.character(grp))
                stop("grp must be either a character vector or a vector of length = ",
                     ncol(object))
              else
                grp <- object[[grp]]

            if(log.it)
              xdat <- log2(object)
            else
              xdat <- object
            

            ans <- stat3(xdat,grp)
            
            
            grp <- factor(grp)
            
            design <- model.matrix(~0+grp)
            colnames(design) <- levels(grp)
            
            ## From LIMMA
            fit <- lm.series(xdat, design = design)
            
            
            ans <- c(ans,list(coefficients=fit$coefficients,
                              cov.coefficients=fit$cov.coefficients,
                              stdev.unscaled=fit$stdev.unscaled))
            
            return(ans)
          }
          )


setMethod("mwt2",
    signature(object="ExpressionSet"),
    function (object, grp, log.it=FALSE) 
          {
            
            if(length(grp) == 1)
              if(!is.character(grp))
                stop("grp must be either a character vector or a vector of length = ",
                     ncol(exprs(object)))
              else
                if(!grp %in% names(pData(object)))
                  stop("grp must be an existing variable in object's phenoData")
                else
                  grp <- object[[grp]]

            if(log.it)
              xdat <- log2(exprs(object))
            else
              xdat <- exprs(object)
            
            
            ans <- stat3(xdat,grp)
            
            grp <- factor(grp)
            
            design <- model.matrix(~0+grp)
            colnames(design) <- levels(grp)
            
            ## From LIMMA
            fit <- lm.series(xdat, design = design)
            
            ans <- c(ans,list(coefficients=fit$coefficients,
                              cov.coefficients=fit$cov.coefficients,stdev.unscaled=fit$stdev.unscaled))
            
            return(ans)
          }
          )



## F statistic for 2 groups
## 24 Feb 2009 - Now uses Levene test to get fFDR
## 7 June 2016 - Add local local FDR

stat2 <- function(xdat,grp,locfdr=FALSE,na.rm=TRUE)
{
    ## basic statistics
    glab = unique(grp)
    n1 = sum(grp==glab[1])
    n2 = sum(grp==glab[2])
    d1 = n1-1
    d2 = n2-1
    m1 = rowMeans(xdat[,grp==glab[1]], na.rm=na.rm)
    m2 = rowMeans(xdat[,grp==glab[2]], na.rm=na.rm)
    
    s2.g1 = rowSums((xdat[,grp==glab[1]]-m1)^2, na.rm=na.rm)/d1
    
    ## We might either have all NA in one group or variance = 0
    ## (e.g. might happen with RMA with small samples)
    ## In this situation we want to remove the gene
    s2.g1[s2.g1 == 0] <- NA
    
    s2.g2 = rowSums((xdat[,grp==glab[2]]-m2)^2, na.rm=na.rm)/d2
    s2.g2[s2.g2 == 0] <- NA
    
    ## If either s2.g1 or s2.g2 are NA this will be NA
    sig2 = (d1*s2.g1 + d2*s2.g2)/(d1+d2)
    fac = 1/n1 + 1/n2
    se2 = (sig2 * fac)
    
    ## F test
    
    lev.test = levene(xdat,grp)
    fFDR = lev.test$FDR
    fStat = lev.test$statistic
    
    
    ## ordinary Welch statistics
    se2.sep = s2.g1/n1 + s2.g2/n2
    df = se2.sep^2/((s2.g1/n1)^2/d1 + (s2.g2/n2)^2/d2)
    
    ## weighted formulas
    df.w  = fFDR*(d1+d2) + (1-fFDR)*df
    se2.w = fFDR*se2 + (1-fFDR)*se2.sep
    ds = est.hyper(z=log(se2.w),D=mean(df.w,na.rm=na.rm),d12=d1+d2)   
    
    ## ....................................... moderated Welch
    se2.com = (ds$d0*ds$s2 + df.w*se2.w)/(ds$d0 + df.w)
    Wm = (m1-m2)/sqrt(se2.com) ## Welch t
    df.com = ds$d0 + df.w      ## df
    

    Wm.pval = pt(-abs(Wm), df= df.com) * 2
    
    ## ................. Compute Global FDR

    Wm.FDR = pval2FDR(Wm.pval) ## Global FDR

    p0_mwt = attr(Wm.FDR,"p0")
    
    ## ................ Compute local FDR

    fdr <- NULL
    
    if(locfdr)
    {
        nr = 50

        Lf.pval = pt(Wm, df= df.com)
        Z = qnorm(Lf.pval)
        
        xbreaks = OCplus:::MidBreaks(Z, nr)
        xmids = brk2mid(xbreaks)
        
        ## ... Smooth fz
        count = hist(Z, xbreaks, plot = FALSE)$counts
        scount = OCplus:::smooth1d(count, err = 1e-04, sv2 = 0.01, verb = FALSE)$fit
        fz = scount/(length(Z) * (xbreaks[2] - xbreaks[1]))
        
        ## ... for f0 let's use Standard Normal density
        f0 <- dnorm(xmids)
        f0fz = f0/fz
        
        ## zlim = 1
        ## p0_OC = 1/max(f0fz[abs(xbreaks) < zlim],na.rm=TRUE)

        sfdr = p0_mwt * f0fz
        fdr = approx(xmids, sfdr, xout = Z, rule = 2)$y
    }

    return(list(MWT= Wm, coefficients=cbind((m1-m2)),pvalue = Wm.pval,
                FDR = Wm.FDR, fdr=fdr,
                df=df.com, s2.wm=se2.com, d0.prior = ds$d0,p0 = p0_mwt,
                s2.prior = ds$s2,lev.stat = fStat, lev.FDR=fFDR))
    
}


## F statistic for 3 or more groups
## Welch statistic + BF statistic
stat3 <- function(xdat,grp,na.rm=TRUE)
{
  glab = unique(grp)
  ngr = length(glab)
  n = mn = s2 = NULL
  
  for (i in 1:ngr)
    { 
      ndx = grp == glab[i]
      ni = sum(ndx)
      mni = rowMeans(xdat[, ndx],na.rm=na.rm)
      x0 = xdat[, ndx] - mni
      s2i = rowSums(x0 * x0, na.rm = na.rm)/(ni-1)

      s2i[s2i == 0] <- NA

      n = c(n, ni)
      mn = cbind(mn, mni)
      s2 = cbind(s2, s2i)
    } 

  ## .............   F statistic
  N = sum(n)
  mmn = rowSums(xdat,na.rm=na.rm)/N
  mn0 = mn - mmn
  num = colSums(t(mn0*mn0)*n,na.rm=na.rm)/(ngr-1)
  den = colSums(t(s2)*(n-1),na.rm=na.rm)/(N-ngr)
  ## F3 = num/den
  df1 = ngr-1

  ## for the weighted sum of squares
  ss.eq = (ngr-1)*den  # pooled sum of squares

  ## ............    Brown-Forsyth
  den = colSums(t(s2)*(1-n/N),na.rm=na.rm)
  ss.uneq = den
  BF = (ngr-1)*num/den
  ci = t(t(s2)*(1-n/N))   # 
  ci = ci/rowSums(ci,na.rm=na.rm)
  df2.BF = 1/colSums(t(ci*ci)/(n-1),na.rm=na.rm)
  
  
  ## weighted SS, df and modified F

  lev.test <- levene(xdat,grp)
  fFDR <- lev.test$FDR
  fStat <- lev.test$statistic
  
  ss.w = fFDR*ss.eq + (1-fFDR)*ss.uneq
  df.w  = fFDR*(N-ngr) + (1-fFDR)*df2.BF
  ds = est.hyper(z=log(ss.w),D=mean(df.w,na.rm=na.rm),d12=N-ngr)   # D = mean(df) 

  df2.Fm = ds$d0 + df.w
  
  ss.m = (ds$d0*ds$s2 + df.w*ss.w)/df2.Fm
  Fm = (ngr-1)*num/ss.m
  
  fpval = 1-pf(Fm, df1=df1, df2 = df2.Fm)
  Fm.FDR = pval2FDR(fpval)
  
  
  return(list("moderated.F" = Fm,df.num=df1, s2.Fm = ss.m/df1,
              df.den = df2.Fm, moderated.F.FDR = Fm.FDR,lev.stat = fStat,
              lev.FDR=fFDR, d0.prior = ds$d0,
              s2.prior = ds$s2))
}





## From contrast.fit in limma
contrasts.mwt <- 
  function (fit, contrasts = NULL, coefficients = NULL) 
{
  ncoef <- NCOL(fit$coefficients)
  if (is.null(contrasts) == is.null(coefficients)) 
    stop("Must specify only one of contrasts or coefficients")
  if (!is.null(contrasts)) {
    rn <- rownames(contrasts)
    cn <- colnames(fit$coefficients)
    if (!is.null(rn) && !is.null(cn) && any(rn != cn)) 
      warning("row names of contrasts don't match col names of coefficients")
  }
  if (!is.null(coefficients)) {
    ncont <- length(coefficients)
    contrasts <- diag(ncoef)
    rownames(contrasts) <- colnames(contrasts) <- colnames(fit$coefficients)
    contrasts <- contrasts[, coefficients, drop = FALSE]
  }
  if (NROW(contrasts) != ncoef) 
    stop("Number of rows of contrast matrix must match number of coefficients")
  fit$contrasts <- contrasts
  cormatrix <- cov2cor(fit$cov.coefficients)
  if (is.null(cormatrix)) {
    warning("no coef correlation matrix found in fit - assuming orthogonal")
    cormatrix <- diag(ncoef)
  }
  r <- nrow(cormatrix)
  if (r < ncoef) {
    if (is.null(fit$pivot)) 
      stop("cor.coef not full rank but pivot column not found in fit")
    est <- fit$pivot[1:r]
    if (any(contrasts[-est, ])) 
      stop("trying to take contrast of non-estimable coefficient")
    contrasts <- contrasts[est, , drop = FALSE]
    fit$coefficients <- fit$coefficients[, est, drop = FALSE]
    fit$stdev.unscaled <- fit$stdev.unscaled[, est, drop = FALSE]
    ncoef <- r
  }


  na.id <- is.na(fit$coefficients)
  fit$stdev.unscaled[na.id] <- 0
  fit$coefficients[na.id] <- 0
  na.id <- na.id %*% contrasts
  na.id[na.id != 0] <- 1
  mode(na.id) <- "logical"
  fit$coefficients <- fit$coefficients %*% contrasts
  fit$coefficients[na.id] <- NA
  
  if (length(cormatrix) < 2) {
    orthog <- TRUE
  }
  else {
    orthog <- all(abs(cormatrix[lower.tri(cormatrix)]) < 
                  1e-14)
  }
  R <- chol(fit$cov.coefficients)
  fit$cov.coefficients <- crossprod(R %*% contrasts)
  fit$pivot <- NULL
  if (orthog)
    {
      fit$stdev.unscaled <- sqrt(fit$stdev.unscaled^2 %*% contrasts^2)
      fit$stdev.unscaled[na.id] <- NA
    }
  else {
    R <- chol(cormatrix)
    ngenes <- NROW(fit$stdev.unscaled)
    ncont <- NCOL(contrasts)
    U <- matrix(1, ngenes, ncont, dimnames = list(rownames(fit$stdev.unscaled), 
                                    colnames(contrasts)))
    o <- array(1, c(1, ncoef))
    for (i in 1:ngenes) {
      RUC <- R %*% limma:::.vecmat(fit$stdev.unscaled[i, ], contrasts)
      U[i, ] <- sqrt(o %*% RUC^2)
    }
    fit$stdev.unscaled <- U
    fit$stdev.unscaled[na.id] <- NA
  }

  fit$t <- fit$coefficients/fit$stdev.unscaled/sqrt(fit$s2.Fm)
    
  t.p.value <- 2 * pt(-abs(fit$t), df = fit$df.den)
  t.fdr <- apply(t.p.value,2,pval2FDR)
  fit$p.value <- t.p.value
  fit$FDR <- t.fdr
  
  fit
}
