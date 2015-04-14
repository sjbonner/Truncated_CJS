fitCJS <- function(W,T,k,control=NULL,hess=FALSE,truth=NULL,debug=FALSE){
    if(debug)
        browser()
    
    ## Compute summary statistics
    R <- summstatCJS(W,k)
    
    ## Compute MLEs
    llhd.opt <- mleCJS(R,T,k,rep(.5,T-1),rep(.5,T-2),control=control,hess=hess)
    
    ## Results
    ## 1) MLEs
    eta.phi.hat <- llhd.opt$par[1:(T-1)]
    phi.hat <- logit$linkinv(eta.phi.hat)
    
    eta.p.hat <- llhd.opt$par[(T-1)+(1:(T-1))]
    p.hat <- logit$linkinv(eta.p.hat)
    
    ## 2) Confidence intervals
    if(!is.null(llhd.opt$hessian)){
        cov <- -solve(llhd.opt$hessian)
        
        eta.phi.sd <- sqrt(diag(cov))[1:(T-1)]
        phi.ci <- logit$linkinv(eta.phi.hat + 2*outer(eta.phi.sd,c(-1,1)))
        
        eta.p.sd <- sqrt(diag(cov))[1:(T-1)]
        p.ci <- logit$linkinv(eta.p.hat + 2*outer(eta.p.sd,c(-1,1)))

        phi.out <- cbind(phi.hat,phi.ci)
        p.out <- cbind(p.hat,p.ci)
    }
    else{
        phi.out <- phi.hat
        p.out <- p.hat
    }

    if(!is.null(truth)){
        phi.out <- cbind(truth$phi,phi.out)
        p.out <- cbind(truth$p,p.out)
    }

    list(phi=phi.out,
         p=p.out,
         opt=llhd.opt)
    
}
    
