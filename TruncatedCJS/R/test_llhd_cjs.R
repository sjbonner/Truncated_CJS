
testdP.dphi <- function(x,t,s,pars,T,k){
    phi <- x
    p <- c(logit$linkinv(pars[T:(2*T-3)]),1)
    #p <- logit$linkinv(pars[T:(2*T-2)])

    PCJS(phi,p,T,k)[t,s]
}

testdP.dp <- function(x,t,s,pars,T,k){
    phi <- logit$linkinv(pars[1:(T-1)])
    p <- c(logit$linkinv(pars[T:(2*T-3)]),1)

    p <- c(x,1)

    P <- PCJS(phi,p,T,k)[t,s]
}
    

testd2P.dphi2 <- function(x,u,v,t,s,pars,T,k){
    phi <- logit$linkinv(pars[1:(T-1)])
    p <- c(logit$linkinv(pars[T:(2*T-3)]),1)

    phi[u] <- x

    P <- PCJS(phi,p,T,k)
    
    dP.dphiCJS(phi,p,P,T,k)[t,s,v]
}

testd2P.dp2 <- function(x,u,v,t,s,pars,T,k){
    phi <- logit$linkinv(pars[1:(T-1)])
    p <- c(logit$linkinv(pars[T:(2*T-3)]),1)

    p[u] <- x

    P <- PCJS(phi,p,T,k)
    
    dP.dpCJS(phi,p,P,T,k)[t,s,v]
}

testd2P.dphidp <- function(p,u,t,s,pars,T,k){
    
    phi <- logit$linkinv(pars[1:(T-1)])
    p <- c(p,1)
    
    P <- PCJS(phi,p,T,k)
    
    tmp <- dP.dphiCJS(phi,p,P,T,k)[t,s,u]

    if(is.na(tmp))
        browser()

    return(tmp)

}

testl.phi <- function(phi,p,R,T,k){
    pars <- c(logit$linkfun(phi),logit$linkfun(p[1:(T-2)]))

    llhdCJS(pars,R,T,k)
}

testl.p <- function(p,phi,R,T,k){

    if(any(phi<0) || any(phi>1))
        return(-Inf)
    
    if(any(p<0) || any(p>1))
        return(-Inf)
    
    pars <- c(logit$linkfun(phi),logit$linkfun(p))

    llhdCJS(pars,R,T,k)
}
