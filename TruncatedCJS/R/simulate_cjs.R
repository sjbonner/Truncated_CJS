simulateCJS <- function(u,T,phi,p,debug=FALSE){

    if(debug)
        browser()

    ## Generate survival matrix
    S <- sapply(1:(T-1),function(t){
        S.tmp <- matrix(NA,u[t],T)
        S.tmp[,t] <- 1
        
        for(s in t:(T-1))
            S.tmp[,s+1] <- 1*(runif(u[t]) < S.tmp[,s] * phi[s])

        S.tmp
    })

    ## Generate capture matrix
    W <- do.call("rbind",sapply(1:(T-1),function(t){
        
        W.tmp <- matrix(0,u[t],T)
        W.tmp[,t] <- 1

        for(s in t:(T-1))
            W.tmp[,s+1] <- 1*(runif(u[t]) < S[[t]][,s+1] * p[s])

        W.tmp
    }))

    return(W)
}
