## Packages
library(parallel)

## Source files
source("R/simulate_cjs.R")
source("R/summary_statistics_cjs.R")
source("R/llhd_cjs.R")
source("R/mle_cjs.R")
source("R/fit_cjs.R")

## Utilities
logit <- make.link("logit")

## Simulation parameters
## nickname <- "Good"
nickname <- "Bad"
## nickname <- "Ugly"

nreps <- 1000
basedir <-paste0("~/Tmp/Big_Data/CJS/",nickname)
outdir <- file.path(basedir,"Output")
datadir <- file.path(basedir,"Data")
parsdir <- file.path(basedir,"Parameters")

if(!file.exists(basedir))
    dir.create(basedir,recursive=TRUE)
if(!file.exists(outdir))
    dir.create(outdir)
if(!file.exists(datadir))
    dir.create(datadir)
if(!file.exists(parsdir))
    dir.create(parsdir)

## Simulation function
runsim <- function(id){
    ## Good model parameters
    if(nickname=="Good"){
        T <- 10    
        phi <- round(rbeta(T-1,3,6),2)
        p <- round(rbeta(T-1,6,2),2)
    }
    if(nickname=="Bad"){
        T <- 10
        phi <- round(rbeta(T-1,6,2),2)
        p <- round(rbeta(T-1,3,6),2)
    }        
    if(nickname=="Ugly"){
        T <- 10
        phi <- rep(.9,T-1) #round(rbeta(T-1,6,2),2)
        p <- c(rep(.1,T-2),.9) #round(rbeta(T-1,3,6),2)
    }        

    truth <- list(phi=phi,p=p)
    
    save(T,phi,p,file=file.path(parsdir,paste0("pars_",id,".Rdata")))
    
    ## Number of unmarked individuals captured on each occasion
    u <- rnbinom(T-1,25,.5)

    ## Simulate data
    W <- simulateCJS(u,T,phi,p)

    ## Store data
    save(W,file=file.path(datadir,paste0("data_",id,".Rdata")))

    ## Fit model
    cat(id,": k=2...\n",sep="")
    fit2 <- fitCJS(W,T,2,control=list(trace=0,maxit=10000),truth=truth)
    save(fit2,file=file.path(outdir,paste0("results_",id,"_",2,".Rdata")))
    
    cat(id,": k=3...\n",sep="")
    fit3 <- fitCJS(W,T,3,control=list(trace=0,maxit=10000),truth=truth)
    save(fit3,file=file.path(outdir,paste0("results_",id,"_",3,".Rdata")))

    cat(id,": k=",T,"...\n",sep="")
    fitT<- fitCJS(W,T,T-1,control=list(trace=0,maxit=10000),truth=truth)
    save(fitT,file=file.path(outdir,paste0("results_",id,"_",T,".Rdata")))

    return("Done")
}

## Run simulations
mclapply(1:nreps,runsim,mc.cores=6)

