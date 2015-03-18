markdataCJS <- function(W,k,file){
    ## Create mark data set

    ## Useful constants
    T <- ncol(W)

    ## Initialize data file
    cat("/*",date(),"*/\n\n",file=file)

    for(i in 1:nrow(W)){
        ## Identify releases
        rels <- which(W[i,]==1)

        if(length(rels)>1){
            for(j in 1:(length(rels)-1)){
                ## Initialize mark data string
                wout <- rep(0,T)
                
                t <- rels[j]
                wout[t] <- 1
                
                wout[(t+1):(rels[j+1]-1)] <- W[i,(t+1):(rels[j+1]-1)]
                wout[rels[j+1]] <- 1
                
                if(rels[j+1]<T)
                    wout[(rels[j+1]+1):T] <- "."

                cat(paste(wout,collapse="")," 1; /*",i,",",t,"*/\n",sep="",file=file,append=TRUE)
            }
        }

        t <- rels[length(rels)]
        if(t<T){
            wout <- rep(0,T)
            wout[t] <- 1
            wout[(t+1):min(t+k,T)] <- rep(0,min(k,T-t))
            if(t+k<T)
                wout[(min(t+k,T)+1):T] <- "."
            cat(paste(wout,collapse="")," 1; /*",i,",",t,"*/\n",sep="",file=file,append=TRUE)
        }
    }
}
