markdataCJSfull <- function(W,file){
    ## Create mark data set for full CJS model

    ## Useful constants
    T <- ncol(W)

    ## Initialize data file
    cat("/*",date(),"*/\n\n",file=file)
    #system(paste("unix2dos",file))

    ## Write capture histories
    for(i in 1:nrow(W))o
        cat(paste(W[i,],collapse="")," 1; /*",i,",",t,"*/\n",sep="",file=file,append=TRUE)
}


