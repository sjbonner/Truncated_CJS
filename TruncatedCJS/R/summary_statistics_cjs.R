summstatCJS <- function(W,k){
    ## Generate summary statistics
    ## k is the number of capture occasions to look ahead

    R <- t(sapply(1:(T-1),function(t){
        ## Identify releases
        rels <- which(W[,t]==1)

        ## Initialze vector of recapture times
        r <- rep(0,k+2)
        r[1] <- length(rels)
        
        ## Identify time of first capture
        for(s in 1:min(k,T-t)){
            recaps <- which(W[rels,t+s]==1)
            r[s+1] <- length(recaps)
            if(length(recaps) > 0)
                rels <- rels[-recaps]
        }
        r[k+2] <- length(rels)

        r
    }))
}    
       
