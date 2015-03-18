## Simulation parameters
nreps <- 1000
basedir <-"~/Tmp/Big_Data/CJS/Bad"
outdir <- file.path(basedir,"Output")
datadir <- file.path(basedir,"Data")
pardir <- file.path(basedir,"Parameters")

## Load output
output <- lapply(1:nreps,function(i){
    file <- file.path(outdir,paste0("results_",i,"_2.Rdata"))
    if(file.exists(file))
        load(file)
    else
        fit2 <- NULL

    file <- file.path(outdir,paste0("results_",i,"_3.Rdata"))
    if(file.exists(file))
        load(file)
    else
        fit3 <- NULL

    file <- file.path(outdir,paste0("results_",i,"_10.Rdata"))
    if(file.exists(file))
        load(file)
    else
        fitT <- NULL

    list(fit2,fit3,fitT)
})

## Extract errors
phi.errors <- sapply(1:nreps,function(i){

    if(is.null(output[[i]][[1]]))
        err2 <- rep(NA,T-1)
    else
        err2 <- output[[i]][[1]]$phi[,2]-output[[i]][[1]]$phi[,1]

    if(is.null(output[[i]][[2]]))
        err3 <- rep(NA,T-1)
    else
        err3 <- output[[i]][[2]]$phi[,2]-output[[i]][[2]]$phi[,1]

    if(is.null(output[[i]][[3]]))
        errT <- rep(NA,T-1)
    else
        errT <- output[[i]][[3]]$phi[,2]-output[[i]][[3]]$phi[,1]

    cbind(err2,err3,errT)
},simplify="array")

p.errors <- sapply(1:nreps,function(i){

    if(is.null(output[[i]][[1]]))
        err2 <- rep(NA,T-1)
    else
        err2 <- output[[i]][[1]]$p[,2]-output[[i]][[1]]$p[,1]

    if(is.null(output[[i]][[2]]))
        err3 <- rep(NA,T-1)
    else
        err3 <- output[[i]][[2]]$p[,2]-output[[i]][[2]]$p[,1]

    if(is.null(output[[i]][[3]]))
        errT <- rep(NA,T-1)
    else
        errT <- output[[i]][[3]]$p[,2]-output[[i]][[3]]$p[,1]

    cbind(err2,err3,errT)
},simplify="array")

## Compare bias
phi.bias <- apply(phi.errors,c(1,2),mean,na.rm=TRUE)
p.bias <- apply(p.errors,c(1,2),mean,na.rm=TRUE)

## 1) Phi
matplot(phi.bias,type="b",ylim=c(-1,1)*max(abs(phi.bias),na.rm=TRUE))
abline(h=0)

## 2) p
matplot(p.bias,type="b",ylim=c(-1,1)*max(abs(p.bias),na.rm=TRUE))
abline(h=0)

## Compare variability
phi.var <- apply(phi.errors,c(1,2),var,na.rm=TRUE)
p.var <- apply(p.errors,c(1,2),var,na.rm=TRUE)

## 1) Phi
matplot(phi.var,type="b",ylim=c(0,max(phi.var,na.rm=TRUE)))
abline(h=0)

## 2) p
matplot(p.var,type="b",ylim=c(0,max(p.var,na.rm=TRUE)))
abline(h=0)

## Compare MSE
phi.mse <- apply(phi.errors^2,c(1,2),mean,na.rm=TRUE)
p.mse <- apply(p.errors^2,c(1,2),mean,na.rm=TRUE)

## 1) Phi
matplot(phi.mse,type="b",ylim=c(0,max(phi.mse,na.rm=TRUE)))
abline(h=0)

## 2) p
matplot(p.mse,type="b",ylim=c(0,max(p.mse,na.rm=TRUE)))
abline(h=0)


