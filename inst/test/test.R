library(mwt)
library(multcomp)

d0=10; s2 = 1   
ng = 100
p0=0.95; 

eqvar=FALSE

n2= n3 = n1 = 3

xdat = MAsim.uneqvar(ng = ng, n1 = n1, n2 = n2, p0 = p0, d01 = 10, s2_01 = 1, v0m = 1, v0var = 1, p0var = 0.7)

grp = colnames(xdat)

a2 <- mwt(xdat,grp)

a3 <- mwt2(xdat,grp)

cMat <- cbind("1-0"=c(-1,1))
rownames(cMat) = c("0","1")


a3.1 <- contrasts.mwt(a3,cMat)

## Brown-Forsythe test and Welch test are equivalent when k=2
all.equal(abs(a3.1$coef),abs(a2$coef),check=FALSE) ## OK
all.equal(a3.1$FDR[,1],a2$FDR) ## OK
all.equal(a3.1$moderated.F,a2$MWT^2) ## OK

#######


xdat3 = MAsim.var6(ng = 1000, n1 = n1, n2 = n2,n3=n3, p0 = p0, d01 = 10, s2_01 = 1, v0m = 1,
  v0var = 1, p0var = 0.7)
grp3 = colnames(xdat3)

a4 <- mwt2(xdat3,grp3)

cMat <- t(contrMat(table(grp3),"Tukey"))
a4.1 <- contrasts.mwt(a4,cMat)

colSums(a4.1$p.value < 0.05)
colSums(a4.1$FDR < 0.05)

