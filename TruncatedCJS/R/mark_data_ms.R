markdataMSind <- function(w,k){
    ## Create capture histories for each individual

    ## Useful constants
    T <- length(w)

    ## Identify releases
    caps <- which(w>0)
    ncaps <- length(caps)

    ## Initialize output matrix
    Wout <- matrix(".",sum(caps<T),T)

    if(ncaps>1){
        for(j in 1:(ncaps-1)){
            ## Initialize mark data string
            t <- caps[j]
            Wout[j,t] <- w[t]

            if(caps[j+1]-1>caps[j])
                Wout[j,(t+1):min(caps[j+1]-1,min(T,t+k))] <- 0

            if(caps[j+1]-caps[j] <= k)
                Wout[j,caps[j+1]] <- w[caps[j+1]]

            ## if(caps[j+1]<T)
            ##     Wout[j,(caps[j+1]+1):T] <- "."

        }
    }

    t <- caps[ncaps]
    if(t<T){
        Wout[ncaps,t] <- w[t]
        Wout[ncaps,(t+1):min(t+k,T)] <- rep(0,min(k,T-t))
        if(t+k<T)
                Wout[ncaps,(min(t+k,T)+1):T] <- "."
    }

    Wout
}

markdataMS <- function(W,k,file=NULL,freq=NULL){

    ## If history frequencies aren't provided then set to 1
    if(is.null(freq))
        freq <- rep(1,nrow(W))

    ## Create truncated capture histories
    Wnew.list <- apply(W,1,markdataMSind,k=k)
    freq.new <- rep(freq,sapply(Wnew.list,nrow))
    Wnew <- do.call("rbind",Wnew.list)

    ## Regroup data
    release <- apply(Wnew,1,function(w) min(which(w!=".")))

    markdata <- do.call("rbind",lapply(1:(T-1),
        function(r){
            tmp <- which(release==r)

            Wunique <- unique(Wnew[tmp,])

            freq.cum <- apply(Wunique,1,function(w1){
                tmp2 <- which(apply(Wnew[tmp,],1,function(w2) all(w1==w2)))
                sum(freq.new[tmp2])
            })

            data.frame(ch=apply(Wunique,1,paste,collapse=""),freq=freq.cum,
                       stringsAsFactors=FALSE)
        }))

    if(!is.null(file)){
        cat("/*",date(),"*/\n\n",file=file)
        write.table(markdata,file,append=TRUE,quote=FALSE,row.names=FALSE,col.names=FALSE)
    }

    return(markdata)

}


