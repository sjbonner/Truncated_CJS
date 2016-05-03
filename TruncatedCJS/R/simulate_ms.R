simulateMS <- function(nstate,u,T,phi,p,debug=FALSE){

    if(debug)
        browser()

    ## Generate complete state matrix
    Wfull <- sapply(1:(T-1),function(t){
        Wfull.tmp <- matrix(NA,sum(u[,t]),T)
        Wfull.tmp[,t] <- rep(1:nstate,u[,t])

        for(s in t:(T-1)){
            for(m in 1:nstate){
                tmp <- which(Wfull.tmp[,s]==m)

                if(length(tmp) > 0)
                    Wfull.tmp[tmp,s+1] <- sample(nstate,length(tmp),replace=TRUE,psi[s,m,])
            }
        }

        Wfull.tmp
    })

    ## Generate survival matrix
    S <- sapply(1:(T-1),function(t){
        S.tmp <- matrix(NA,sum(u[,t]),T)
        S.tmp[,t] <- 1

        for(s in t:(T-1))
            S.tmp[,s+1] <- 1*(runif(sum(u[,t])) < S.tmp[,s] * phi[Wfull[[t]][,s],s])

        S.tmp
    })

    ## Generate capture matrix
    W <- do.call("rbind",sapply(1:(T-1),function(t){

        W.tmp <- matrix(0,sum(u[,t]),T)
        W.tmp[,t] <-Wfull[[t]][,t]

        
        for(s in t:(T-1)){
            p.tmp <- S[[t]][,s+1] * p[Wfull[[t]][,s+1],s]
            W.tmp[,s+1] <- Wfull[[t]][,s+1]*(runif(sum(u[,t])) <  p.tmp)
        }

        W.tmp
    }))

    return(W)
}
