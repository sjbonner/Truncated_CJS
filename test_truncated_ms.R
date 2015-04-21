## Clear the workspace of cruft
rm(list=ls())

## Load devtools package
library(devtools)

## Load package from library
## library(TruncatedCJS)

## Reload files from package directory
## devtools::reload(packdir)

## Load all package files
packdir <- "TruncatedCJS"
devtools::load_all(packdir)

## Set seed
set.seed(777)

## Constant time model
nstate <- 5

usize <- c(25,50,100,100,50)
T <- 10

phi <- matrix(rep(round(rbeta(nstate,12,4),2),T-1),nrow=nstate)
p <- matrix(rep(round(rbeta(nstate,12,4),2),T-1),nrow=nstate)

psi <- array(NA,dim=c(T-1,nstate,nstate))

psi.tmp <- c(.3,rep(.7/4,4))

for(t in 1:(T-1)){
    for(m1 in 1:nstate){
        for(m2 in 1:nstate){
            psi[t,m1,m2] <- psi.tmp[abs(m1-m2)+1]
        }
    }
}

truth <- list(phi=phi,p=p,psi=psi)

## Number of unmarked individuals captured in each strata on each occasion
u <- matrix(rnbinom(nstate*(T-1),usize,.5),nrow=nstate)

## Simulate data
W <- simulateMS(nstate,u,T,phi,p,debug=FALSE)

## Generate MARK data
k <- 3 #T-1
msdata <- markdataMS(W,k,paste0("MARK/mark_ms_test_",k,".INP"))

## Fit multistate model in RMark
msdata.processed <- process.data(msdata,model="Multistrata")

## Stratum effects only
msdata.ddl1 <- make.design.data(msdata.processed,
                                          parameter=list(Psi=list(pim.type="constant"),
                                              S=list(pim.type="constant"),
                                              p=list(pim.type="constant")))

msdata.ddl1$Psi$d <- as.factor(abs(as.numeric(msdata.ddl1$Psi$stratum)-
                                   as.numeric(msdata.ddl1$Psi$tostratum)))

## Difference based model of transitions
## formd <- list(formula= ~ -1+stratum:d)

## time1 <- system.time(ms.model1 <- mark(msdata.processed,msdata.ddl1,threads=3,
##                                        model.parameters=list(Psi=formd)))

## Unconstrained model
time1 <- system.time(ms.model1 <- mark(msdata.processed,msdata.ddl1,threads=3))

## Plot results

## 1) Phi
par(mfrow=c(1,1))
index <- 1:nstate

plot(1:nstate,phi[,1],pch=16,col="blue",ylim=c(0,1))
points(1:nstate,ms.model1$results$real[index,"estimate"],pch=16)
for(i in 1:nstate)
    lines(rep(i,2),ms.model1$results$real[index[i],c("lcl","ucl")])

## 2) p
par(mfrow=c(1,1))
index <- nstate + 1:nstate

plot(1:nstate,p[,1],pch=16,col="blue",ylim=c(0,1))
points(1:nstate,ms.model1$results$real[index,"estimate"],pch=16)
for(i in 1:nstate)
    lines(rep(i,2),ms.model1$results$real[index[i],c("lcl","ucl")])

## 3) psi
par(mfrow=c(ceiling(nstate/2),2))

for(m in 1:nstate){
    index <- 2 * nstate + (m-1) * (nstate-1) + 1:(nstate-1)

    plot((1:nstate),psi[1,m,],pch=16,col="blue",ylim=c(0,1),xlim=c(1,nstate))
    points((1:nstate)[-m],ms.model1$results$real[index,"estimate"],pch=16)
    points(m,1-sum(ms.model1$results$real[index,"estimate"]),pch=16)
    for(i in 1:(nstate-1)){
        t <- (1:nstate)[-m][i]
        lines(rep(t,2),ms.model1$results$real[index[i],c("lcl","ucl")])
    }
}

