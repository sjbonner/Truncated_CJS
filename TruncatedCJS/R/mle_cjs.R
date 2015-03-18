mleCJS <- function(R,T,k,phi,p,control=NULL,hess=FALSE){
    ## Compute initial values on logistic scale
    inits <- pmax(pmin(c(logit$linkfun(phi),
                         logit$linkfun(p)),
                       logit$linkfun(.999)),
                  logit$linkfun(.001))
    
    ## Maximize likelihood and compute Hessian
    optim(inits,llhdCJS,llhdCJSgr,R=R,T=T,k=k,
          method="BFGS",debug=FALSE,
          control=c(list(fnscale=-1,trace=10),control),hessian=hess)
}
