PCJS <- function(phi,p,T,k){
    ## Define multinomial cell probabilities
    t(sapply(1:(T-1),function(t){
        l <- min(k,T-t)-1

        if(l>0)
            tmp <- cumprod(phi[t:(t+l)]) * cumprod(c(1,1-p[t:(t+l-1)])) * p[t:(t+l)]
        else
            tmp <- phi[t]*p[t]

        c(tmp,rep(0,max(0,k-l-1)),1-sum(tmp,na.rm=TRUE))
    }))
}

dP.dphiCJS <- function(phi,p,P,T,k){
    ## Differentiate cell probabilities wrt survival parameters

    out <- array(0,c(T-1,k,T-1))

    for(t in 1:(T-1)){                  # Released at t
        for(s in 1:min(k,T-t)){         # Recaptured at s
            for(r in t:(t+s-1))
                out[t,s,r] <- P[t,s]/phi[r]
        }
    }

    out
}

d2P.dphi2CJS <- function(phi,p,P,T,k){
    ## Differentiate cell probabilities wrt survival parameters

    out <- array(0,c(T-1,k,T-1,T-1))

    for(t in 1:(T-1)){                  # Released at t
        for(s in 1:min(k,T-t)){         # Recaptured at s
            for(u in t:(t+s-1)){
                for(v in t:(t+s-1)){
                    if(u!=v)
                        out[t,s,u,v] <- P[t,s]/phi[u]/phi[v]
                }
            }
        }
    }

    out
}

dP.dpCJS <- function(phi,p,P,T,k){
    ## Differentiate cell probabilities wrt capture parameters

    out <- array(0,c(T-1,k,T-2))

    for(t in 1:(T-1)){                  # Released at t
        for(s in 1:min(k,T-t)){         # Recaptured at s
            ## Missed captures
            if(s>1){
                for(r in t:(t+s-2)){
                    out[t,s,r] <- -P[t,s]/(1-p[r])
                }
            }

            ## Recapture
            if((t+s-1)<(T-1)){
                out[t,s,t+s-1] <- P[t,s]/p[t+s-1]
            }
        }
    }

    out
}

d2P.dp2CJS <- function(phi,p,P,T,k){
    ## Differentiate cell probabilities wrt capture parameters

    out <- array(0,c(T-1,k,T-2,T-2))

    for(t in 1:(T-1)){                  # Released at t
        if(min(k,T-t)>1){
            for(s in 2:min(k,T-t)){         # Recaptured at s
                for(u in t:(t+s-2)){
                    for(v in t:(t+s-2)){
                        if(u!=v)
                            out[t,s,u,v] <- P[t,s]/(1-p[u])/(1-p[v])
                    }
                }

                if((t+s-1)<(T-1)){
                    for(u in t:(t+s-2))
                        out[t,s,u,t+s-1] <- -P[t,s]/(1-p[u])/p[t+s-1]

                    for(v in t:(t+s-2))
                        out[t,s,t+s-1,v] <- -P[t,s]/p[t+s-1]/(1-p[v])
                }
            }
        }
    }

    out
}

d2P.dphidpCJS <- function(phi,p,P,T,k){
    ## Differentiate cell probabilities wrt survival parameters

    out <- array(0,c(T-1,k,T-1,T-2))

    for(t in 1:(T-1)){                  # Released at t
        for(s in 1:min(k,T-t)){         # Recaptured at s
            for(u in t:(t+s-1)){        # Index on phi
                for(v in t:max(t,(t+s-2))){    # Index on p
                    if(v<T-1){
                        out[t,s,u,v] <- -P[t,s]/phi[u]/(1-p[v])
                    }
                }
                if((t+s-1)<(T-1))
                    out[t,s,u,t+s-1] <- P[t,s]/phi[u]/p[t+s-1]
            }
        }
    }

    out
}

llhdCJS <- function(pars,R,T,k,debug=FALSE){
    if(debug)
        browser()

    ## Unpack parameters
    phi <- logit$linkinv(pars[1:(T-1)])

    if(length(pars)==2*T-3)
        p <- c(logit$linkinv(pars[T:(2*T-3)]),1)
    else
        p <- logit$linkinv(pars[T:(2*T-2)])

    ## Define multinomial probabilities
    P <- PCJS(phi,p,T,k)

    ## Compute likelihood
    sum(sapply(1:(T-1),function(t){
        dmultinom(R[t,-1],R[t,1],P[t,],log=TRUE)
    }))
}

llhdCJSgr <- function(pars,R,T,k,debug=FALSE,logitscale=TRUE){
    if(debug)
        browser()

    ## Unpack parameters
    phi <- logit$linkinv(pars[1:(T-1)])

    if(length(pars)==2*T-3)
        p <- c(logit$linkinv(pars[T:(2*T-3)]),1)
    else
        p <- logit$linkinv(pars[T:(2*T-2)])

    ## Define multinomial probabilities
    P <- PCJS(phi,p,T,k)

    ## Compute derivative of each cell probability wrt phi
    dP.dphi <- dP.dphiCJS(phi,p,P,T,k)

    ## Compute derivative of each cell probability wrt phi
    dP.dp <- dP.dpCJS(phi,p,P,T,k)

    ## Compute derivative of likelihood wrt each element of P
    dl.dP <- R[,2:(k+1)]/P[,1:k] - R[,k+2]/P[,k+1]

    ## Compute gradient
    dl.dphi <- sapply(1:(T-1),function(r){
        sum(dl.dP * dP.dphi[,,r],na.rm=TRUE)
    })

    dl.dp <- sapply(1:(T-2),function(r){
        sum(dl.dP * dP.dp[,,r],na.rm=TRUE)
    })

    if(!logitscale)
        return(c(dl.dphi,dl.dp))

    ## Transform to logit scale
    else
        return(c(dl.dphi,dl.dp) * c(phi,p[-(T-1)])*(1-c(phi,p[-(T-1)])))
}

llhdCJShess <- function(pars,R,T,k,debug=FALSE){
    if(debug)
        browser()

    ## Unpack parameters
    phi <- logit$linkinv(pars[1:(T-1)])

    if(length(pars)==2*T-3)
        p <- c(logit$linkinv(pars[T:(2*T-3)]),1)
    else
        p <- logit$linkinv(pars[T:(2*T-2)])

    ## Define multinomial probabilities
    P <- PCJS(phi,p,T,k)

    ## Compute derivative of each cell probability wrt phi
    dP.dphi <- dP.dphiCJS(phi,p,P,T,k)

    ## Compute derivative of each cell probability wrt p
    dP.dp <- dP.dpCJS(phi,p,P,T,k)

    ## Compute second derivatives of P wrt phi and p
    d2P.dphi2 <- d2P.dphi2CJS(phi,p,P,T,k)
    d2P.dp2 <- d2P.dp2CJS(phi,p,P,T,k)
    d2P.dphidp <- d2P.dphidpCJS(phi,p,P,T,k)

    ## Compute derivative of likelihood wrt each element of P
    dl.dP <- R[,2:(k+1)]/P[,1:k] - R[,k+2]/P[,k+1]

    ## Compute second derivative of likelihood wrt each element of P
    d2l.dP2 <- -R[,2:(k+1)]/P[,1:k]^2 - R[,k+2]/P[,k+1]^2

    ## Compute Hessian
    d2l.dphi2 <- sapply(1:(T-1),function(u){
        sapply(1:(T-1),function(v){
            - sum(R[,2:(k+1)]/P[,1:k]^2 * dP.dphi[,,u] * dP.dphi[,,v],na.rm=TRUE) -
                sum(R[,k+2]/P[,k+1]^2 * apply(dP.dphi[,,u],1,sum) * apply(dP.dphi[,,v],1,sum), na.rm=TRUE) +
                    sum(dl.dP * d2P.dphi2[,,u,v],na.rm=TRUE)
        })
    })

    d2l.dp2 <- sapply(1:(T-2),function(u){
        sapply(1:(T-2),function(v){
            - sum(R[,2:(k+1)]/P[,1:k]^2 * dP.dp[,,u] * dP.dp[,,v],na.rm=TRUE) -
                sum(R[,k+2]/P[,k+1]^2 * apply(dP.dp[,,u],1,sum) * apply(dP.dp[,,v],1,sum), na.rm=TRUE) +
                    sum(dl.dP * d2P.dp2[,,u,v],na.rm=TRUE)
        })
    })

    d2l.dphidp <- t(sapply(1:(T-1),function(u){
        sapply(1:(T-2),function(v){
            - sum(R[,2:(k+1)]/P[,1:k]^2 * dP.dphi[,,u] * dP.dp[,,v],na.rm=TRUE) -
                sum(R[,k+2]/P[,k+1]^2 * apply(dP.dphi[,,u],1,sum) * apply(dP.dp[,,v],1,sum), na.rm=TRUE) +
                    sum(dl.dP * d2P.dphidp[,,u,v],na.rm=TRUE)
        })
    }))

    hess.tmp <- rbind(cbind(d2l.dphi2,d2l.dphidp),
                      cbind(t(d2l.dphidp),d2l.dp2))

    ## Transform to logit scale
    tmp1 <- c(phi,p[-(T-1)])*(1-c(phi,p[-(T-1)]))
    tmp2 <- tmp1*(1-2*c(phi,p[-(T-1)]))

    hess <- hess.tmp * outer(tmp1,tmp1) +
        diag(llhdCJSgr(pars,R,T,k,logitscale=FALSE) * tmp2)
}

llhdCJScombined <- function(pars,R,T,k,debug=FALSE,logitscale=TRUE,llhdscale=1,gradient=TRUE,hessian=TRUE){
    if(debug)
        browser()

    cat(pars,"\n",file="pars_test.R",append=TRUE)

    ## Unpack parameters
    phi <- logit$linkinv(pars[1:(T-1)])

    if(length(pars)==2*T-3)
        p <- c(logit$linkinv(pars[T:(2*T-3)]),1)
    else
        p <- logit$linkinv(pars[T:(2*T-2)])

    ##### Compute likelihood #####

    ## Define multinomial probabilities
    P <- PCJS(phi,p,T,k)

    llhd <- llhdscale*sum(sapply(1:(T-1),function(t){
        dmultinom(R[t,-1],R[t,1],P[t,],log=TRUE)
    }))
    
    ## if(llhd>3000){
    ##     cat(pars,"\n")
    ##     stop(llhd,"\n")
    ## }
    
    ##### Compute gradient #####
    if(gradient || hessian){
        ## Compute derivative of each cell probability wrt phi
        dP.dphi <- dP.dphiCJS(phi,p,P,T,k)
        
        ## Compute derivative of each cell probability wrt phi
        dP.dp <- dP.dpCJS(phi,p,P,T,k)
        
        ## Compute derivative of likelihood wrt each element of P
        dl.dP <- R[,2:(k+1)]/P[,1:k] - R[,k+2]/P[,k+1]
        
        ## Compute gradient
        dl.dphi <- sapply(1:(T-1),function(r){
            sum(dl.dP * dP.dphi[,,r],na.rm=TRUE)
        })
        
        dl.dp <- sapply(1:(T-2),function(r){
            sum(dl.dP * dP.dp[,,r],na.rm=TRUE)
        })

        gr <- c(dl.dphi,dl.dp)
            
        if(!logitscale)
            attr(llhd,"gradient") <- llhdscale * gr
        
        else
            ## Transform to logit scale
            attr(llhd,"gradient") <- llhdscale * gr * c(phi,p[-(T-1)])*(1-c(phi,p[-(T-1)]))
    }
    
    ##### Compute Hessian #####
    ## if(hessian){
    ##     ## Compute second derivatives of P wrt phi and p
    ##     d2P.dphi2 <- d2P.dphi2CJS(phi,p,P,T,k)
    ##     d2P.dp2 <- d2P.dp2CJS(phi,p,P,T,k)
    ##     d2P.dphidp <- d2P.dphidpCJS(phi,p,P,T,k)
        
    ##     ## Compute second derivative of likelihood wrt each element of P
    ##     d2l.dP2 <- -R[,2:(k+1)]/P[,1:k]^2 - R[,k+2]/P[,k+1]^2
        
    ##     ## Compute Hessian
    ##     d2l.dphi2 <- sapply(1:(T-1),function(u){
    ##         sapply(1:(T-1),function(v){
    ##             - sum(R[,2:(k+1)]/P[,1:k]^2 * dP.dphi[,,u] * dP.dphi[,,v],na.rm=TRUE) -
    ##                 sum(R[,k+2]/P[,k+1]^2 * apply(dP.dphi[,,u],1,sum) * apply(dP.dphi[,,v],1,sum), na.rm=TRUE) +
    ##                     sum(dl.dP * d2P.dphi2[,,u,v],na.rm=TRUE)
    ##         })
    ##     })
        
    ##     d2l.dp2 <- sapply(1:(T-2),function(u){
    ##         sapply(1:(T-2),function(v){
    ##             - sum(R[,2:(k+1)]/P[,1:k]^2 * dP.dp[,,u] * dP.dp[,,v],na.rm=TRUE) -
    ##                 sum(R[,k+2]/P[,k+1]^2 * apply(dP.dp[,,u],1,sum) * apply(dP.dp[,,v],1,sum), na.rm=TRUE) +
    ##                     sum(dl.dP * d2P.dp2[,,u,v],na.rm=TRUE)
    ##         })
    ##     })
        
    ##     d2l.dphidp <- t(sapply(1:(T-1),function(u){
    ##         sapply(1:(T-2),function(v){
    ##             - sum(R[,2:(k+1)]/P[,1:k]^2 * dP.dphi[,,u] * dP.dp[,,v],na.rm=TRUE) -
    ##                 sum(R[,k+2]/P[,k+1]^2 * apply(dP.dphi[,,u],1,sum) * apply(dP.dp[,,v],1,sum), na.rm=TRUE) +
    ##                     sum(dl.dP * d2P.dphidp[,,u,v],na.rm=TRUE)
    ##         })
    ##     }))

    ##     hess <- rbind(cbind(d2l.dphi2,d2l.dphidp),
    ##                   cbind(t(d2l.dphidp),d2l.dp2))
        
    ##     ## Transform to logit scale
    ##     if(logitscale){
    ##         tmp1 <- c(phi,p[-(T-1)])*(1-c(phi,p[-(T-1)]))
    ##         tmp2 <- tmp1*(1-2*c(phi,p[-(T-1)]))
        
    ##         hess <- hess * outer(tmp1,tmp1) + diag(gr * tmp2)
    ##     }
        
    ##     attr(llhd,"hessian") <- llhdscale * hess

    ##     dput(gr,"gradient_test.R")
    ##     dput(hess,"hess_test.R")
    ## }
    
    return(llhd)
}

