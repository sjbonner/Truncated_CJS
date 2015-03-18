## Load devtools package
library(devtools)

## Load package from library
library(TruncatedCJS)

## Reload files from package directory
packdir <- "TruncatedCJS"
devtools::reload(packdir)

## Set seed
set.seed(777)

## Good model parameters
T <- 10
phi <- round(rbeta(T-1,3,6),2)
p <- round(rbeta(T-1,6,2),2)

## Bad model parameters
## T <- 10
## phi <- rep(.9,T-1)
## p <- c(rep(.1,T-2),.9)

truth <- list(phi=phi,p=p)

## Utilities
logit <- make.link("logit")

## Number of unmarked individuals captured on each occasion
u <- rnbinom(T-1,100,.5)

## Simulate data
W <- simulateCJS(u,T,phi,p)

## Fit model
fit2 <- fitCJS(W,T,2,control=list(trace=1,maxit=10000),hess=TRUE,truth=truth,debug=FALSE)
hess2 <- llhdCJShess(fit2$opt$par,summstatCJS(W,2),T,2)
max(abs(fit2$opt$hess-hess2))

fit3 <- fitCJS(W,T,3,control=list(trace=1,maxit=10000),hess=TRUE,truth=truth)
hess3 <- llhdCJShess(fit3$opt$par,summstatCJS(W,3),T,3)
max(abs(fit3$opt$hess-hess3))

fitT<- fitCJS(W,T,T,control=list(trace=1,maxit=10000),hess=TRUE,truth=truth)
hessT <- llhdCJShess(fitT$opt$par,summstatCJS(W,T),T,T)
max(abs(fitT$opt$hess-hessT))

## Compare results
col <- grey((2:0)/3)

## Phi
jit <- .2
plot(NA,NA,ylim=c(0,1),xlim=c(0,T))

points(1:(T-1)-jit,fit2$phi[,2],
     pch=16,col=col[1])

points(1:(T-1),fit3$phi[,2],
     pch=16,col=col[2])

points(1:(T-1)+jit,fitT$phi[,2],
     pch=16,col=col[3])

for(t in 1:(T-1)){
    lines(rep(t,2)-jit,fit2$phi[t,3:4],col=col[1])
    lines(rep(t,2),fit3$phi[t,3:4],col=col[2])
    lines(rep(t,2)+jit,fitT$phi[t,3:4],col=col[3])
}

lines(1:(T-1),phi,type="b",col="red",lty=2)

## P
jit <- .2
plot(NA,NA,ylim=c(0,1),xlim=c(0,T))

points(1:(T-1)-jit,fit2$p[,2],
     pch=16,col=col[1])

points(1:(T-1),fit3$p[,2],
     pch=16,col=col[2])

points(1:(T-1)+jit,fitT$p[,2],
     pch=16,col=col[3])

for(t in 1:(T-1)){
    lines(rep(t,2)-jit,fit2$p[t,3:4],col=col[1])
    lines(rep(t,2),fit3$p[t,3:4],col=col[2])
    lines(rep(t,2)+jit,fitT$p[t,3:4],col=col[3])
}

lines(1:(T-1),p,type="b",col="red",lty=2)

## Compare Hessian as a function of k
summary(as.vector(hess3-hess2)^2)
summary(as.vector(hessT-hess3)^2)

summary(diag(solve(hess3)-solve(hess2))^2)
summary(diag(solve(hessT)-solve(hess3))^2)
summary(diag(solve(hess3)-solve(hess2))^2)

## Compare variance estimates
round(-cbind(diag(solve(hess2)),diag(solve(hess3)),diag(solve(hessT))),2)

##### Test Gradient #####
library(numDeriv)

k <- 3
R <- summstatCJS(W,k)

pars <- c(logit$linkfun(phi),logit$linkfun(p[-(T-1)]))

grad1 <- llhdCJSgr(pars,R,T,k)
grad2 <- grad(llhdCJS,pars,R=R,T=T,k=k)

round(tmp <- grad1-grad2,6)
          
max(abs(tmp))

##### Test Hessian #####
k <- 3
R <- summstatCJS(W,k)

pars <- c(logit$linkfun(phi),logit$linkfun(p[-(T-1)]))

system.time(hess1 <- llhdCJShess(pars,R,T,k))
system.time(hess2 <- hessian(llhdCJS,pars,R=R,T=T,k=k))

round(tmp <- hess1-hess2,3)

max(abs(tmp))

round(hess1[1:(T-1),1:(T-1)]-hess2[1:(T-1),1:(T-1)],3)
round(hess1[T:(2*T-3),T:(2*T-3)]-hess2[T:(2*T-3),T:(2*T-3)],3)
round(hess1[1:(T-1),T:(2*T-3)]-hess2[1:(T-1),T:(2*T-3)],3)

##### Run in MARK #####
markdataCJSfull(W,"MARK/mark_sim_test_full.INP")
markdataCJS(W,T-1,"MARK/mark_sim_test_T-1.INP")
markdataCJS(W,2,"MARK/mark_sim_test_2.INP")
markdataCJS(W,3,"MARK/mark_sim_test_3.INP")

##### Test Steps in Gradient Calculation #####

## 1) dP.dphi
phi <- truth$phi; phi[T-1] <- truth$phi[T-t]*truth$p[T-1]
p <- truth$p; p[T-1] <- 1
pars <- c(logit$linkfun(phi),logit$linkfun(p))

P <- PCJS(phi,p,T,k)
dP.dphi <- dP.dphiCJS(phi,p,P,T,k)

for(t in 1:(T-1)){
    for(s in 1:k){

        tmp1 <- grad(testdP.dphi,phi[1:(T-1)],t=t,s=s,pars=pars,T=T,k=k)

        d <- max(abs(tmp1-dP.dphi[t,s,]))
        
        if(d > 1e-6)
            cat(t,s,":",d,"\n")
    }
}

## 2) dP.dp
P <- PCJS(phi,p,T,k)
dP.dp <- dP.dpCJS(phi,p,P,T,k)
pars <- c(logit$linkfun(phi),logit$linkfun(p))

for(t in 1:(T-1)){
    for(s in 1:k){

        tmp1 <- grad(testdP.dp,p[1:(T-2)],t=t,s=s,pars=pars,T=T,k=k)

        d <- max(abs(tmp1-dP.dp[t,s,]))
        
        if(d > 1e-6)
            cat(t,s,":",d,"\n")
    }
}

## 3) dl.dphi
phi <- truth$phi; phi[T-1] <- truth$phi[T-t]*truth$p[T-1]
p <- truth$p; p[T-1] <- 1 #pmin(pmax(truth$p,.85),.15)
pars <- c(logit$linkfun(phi),logit$linkfun(p))

P <- PCJS(phi,p,T,k)
dP.dphi <- dP.dphiCJS(phi,p,P,T,k)

## Compute derivative of likelihood wrt each element of P
dl.dP <- R[,2:(k+1)]/P[,1:k] - R[,k+2]/P[,k+1]

dl.dphi <- sapply(1:(T-1),function(r){
    sum(dl.dP * dP.dphi[,,r],na.rm=TRUE)
})

tmp <- grad(testl.phi,phi,p=p,R=R,T=T,k=k)

round(dl.dphi-tmp,2)

## 4) dl.dp
phi <- truth$phi; phi[T-1] <- truth$phi[T-t]*truth$p[T-1]
p <- truth$p; p[T-1] <- 1
pars <- c(logit$linkfun(phi),logit$linkfun(p))

P <- PCJS(phi,p,T,k)
dP.dp <- dP.dpCJS(phi,p,P,T,k)

## Compute derivative of likelihood wrt each element of P
dl.dP <- R[,2:(k+1)]/P[,1:k] - R[,k+2]/P[,k+1]

dl.dp <- sapply(1:(T-2),function(r){
    sum(dl.dP * dP.dp[,,r],na.rm=TRUE)
})

tmp <- grad(testl.p,p[1:(T-2)],phi=phi,R=R,T=T,k=k)

round(dl.dp-tmp,2)

##### Test Steps in Hessian Computation #####

## 1) d2P.dphi2
phi <- truth$phi; phi[T-1] <- truth$phi[T-t]*truth$p[T-1]
p <- truth$p; p[T-1] <- 1
pars <- c(logit$linkfun(phi),logit$linkfun(p))

P <- PCJS(phi,p,T,k)
d2P.dphi2 <- d2P.dphi2CJS(phi,p,P,T,k)

for(u in 1:(T-1)){
    for(v in 1:(T-1)){
        tmp1 <- t(sapply(1:(T-1),function(t){
            sapply(1:k,function(s){
                grad(testd2P.dphi2,logit$linkinv(pars[u]),t=t,s=s,u=u,v=v,pars=pars,T=T,k=k)
            })
        }))

        d <- max(abs(tmp1-d2P.dphi2[,,u,v]))

        if(d > 1e-6)
            cat(u,v,":",d,"\n")
    }
}

## 2) d2P.dp2
phi <- truth$phi; phi[T-1] <- truth$phi[T-t]*truth$p[T-1]
p <- truth$p; p[T-1] <- 1
pars <- c(logit$linkfun(phi),logit$linkfun(p))

P <- PCJS(phi,p,T,k)
d2P.dp2 <- d2P.dp2CJS(phi,p,P,T,k)

for(u in 1:(T-2)){
    for(v in u:(T-2)){
        tmp1 <- t(sapply(1:(T-1),function(t){
            sapply(1:k,function(s){
                grad(testd2P.dp2,p[u],t=t,s=s,u=u,v=v,pars=pars,T=T,k=k)
            })
        }))

        d <- max(abs(tmp1-d2P.dp2[,,u,v]))

        if(d > 1e-6)
            cat(u,v,":",d,"\n")
    }
}

## 3) d2P.dphidp
phi <- truth$phi; phi[T-1] <- truth$phi[T-t]*truth$p[T-1]
p <- truth$p; p[T-1] <- 1
pars <- c(logit$linkfun(phi),logit$linkfun(p))

P <- PCJS(phi,p,T,k)
d2P.dphidp <- d2P.dphidpCJS(phi,p,P,T,k)

for(t in 1:(T-1)){
    for(s in 1:k){
        tmp1 <- t(sapply(1:(T-1),function(u){
            grad(testd2P.dphidp,p[1:8],t=t,s=s,u=u,pars=pars,T=T,k=k)
        }))

        d <- max(abs(tmp1-d2P.dphidp[t,s,,]))

        if(d > 1e-6)
            cat(t,s,":",d,"\n")
    }
}

## 4) d2l.dphi2
phi <- truth$phi; phi[T-1] <- truth$phi[T-t]*truth$p[T-1]
p <- truth$p; p[T-1] <- 1 #pmin(pmax(truth$p,.85),.15)
pars <- c(logit$linkfun(phi),logit$linkfun(p))

P <- PCJS(phi,p,T,k)
dP.dphi <- dP.dphiCJS(phi,p,P,T,k)
d2P.dphi2 <- d2P.dphi2CJS(phi,p,P,T,k)

## Compute derivative of likelihood wrt each element of P
dl.dP <- R[,2:(k+1)]/P[,1:k] - R[,k+2]/P[,k+1]

## Compute second derivative of likelihood wrt each element of P
#d2l.dP2 <- -R[,2:(k+1)]/P[,1:k]^2 - R[,k+2]/P[,k+1]^

d2l.dphi2 <- sapply(1:(T-1),function(u){
    sapply(1:(T-1),function(v){
        - sum(R[,2:(k+1)]/P[,1:k]^2 * dP.dphi[,,u] * dP.dphi[,,v],na.rm=TRUE) -
            sum(R[,k+2]/P[,k+1]^2 * apply(dP.dphi[,,u],1,sum) * apply(dP.dphi[,,v],1,sum), na.rm=TRUE) + 
                sum(dl.dP * d2P.dphi2[,,u,v],na.rm=TRUE)
    })
})

tmp <- hessian(testl.phi,phi,p=p,R=R,T=T,k=k)

round(d2l.dphi2-tmp,2)

## 5) d2l.dp2
phi <- truth$phi; phi[T-1] <- truth$phi[T-t]*truth$p[T-1]
p <- pmin(pmax(truth$p,.85),.15); p[T-1] <- 1 # Keep values of p away from boundary
pars <- c(logit$linkfun(phi),logit$linkfun(p))

P <- PCJS(phi,p,T,k)
dP.dp <- dP.dpCJS(phi,p,P,T,k)
d2P.dp2 <- d2P.dp2CJS(phi,p,P,T,k)

## Compute derivative of likelihood wrt each element of P
dl.dP <- R[,2:(k+1)]/P[,1:k] - R[,k+2]/P[,k+1]

## Compute second derivative of likelihood wrt each element of P
#d2l.dP2 <- -R[,2:(k+1)]/P[,1:k]^2 - R[,k+2]/P[,k+1]^

d2l.dp2 <- sapply(1:(T-2),function(u){
    sapply(1:(T-2),function(v){
        - sum(R[,2:(k+1)]/P[,1:k]^2 * dP.dp[,,u] * dP.dp[,,v],na.rm=TRUE) -
            sum(R[,k+2]/P[,k+1]^2 * apply(dP.dp[,,u],1,sum) * apply(dP.dp[,,v],1,sum), na.rm=TRUE) + 
                sum(dl.dP * d2P.dp2[,,u,v],na.rm=TRUE)
    })
})

(tmp <- hessian(testl.p,p[1:8],phi=phi,R=R,T=T,k=k))

round(d2l.dp2-tmp,2)

## 6) d2l.dphidp
phi <- truth$phi; phi[T-1] <- truth$phi[T-t]*truth$p[T-1]
p <- pmin(pmax(truth$p,.85),.15); p[T-1] <- 1 # Keep values of p away from boundary
pars <- c(logit$linkfun(phi),logit$linkfun(p))

P <- PCJS(phi,p,T,k)
dP.dphi <- dP.dphiCJS(phi,p,P,T,k)
dP.dp <- dP.dpCJS(phi,p,P,T,k)
d2P.dphidp <- d2P.dphidpCJS(phi,p,P,T,k)

## Compute derivative of likelihood wrt each element of P
dl.dP <- R[,2:(k+1)]/P[,1:k] - R[,k+2]/P[,k+1]

## Compute second derivative of likelihood wrt each element of P
#d2l.dP2 <- -R[,2:(k+1)]/P[,1:k]^2 - R[,k+2]/P[,k+1]^

d2l.dphidp <- sapply(1:(T-1),function(u){
    sapply(1:(T-2),function(v){
        - sum(R[,2:(k+1)]/P[,1:k]^2 * dP.dphi[,,u] * dP.dp[,,v],na.rm=TRUE) -
            sum(R[,k+2]/P[,k+1]^2 * apply(dP.dphi[,,u],1,sum) * apply(dP.dp[,,v],1,sum), na.rm=TRUE) + 
                sum(dl.dP * d2P.dphidp[,,u,v],na.rm=TRUE)
    })
})

#??(tmp <- hessian(testl.p,p[1:8],phi=phi,R=R,T=T,k=k))

round(d2l.dp2-tmp,2)

