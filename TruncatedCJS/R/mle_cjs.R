mleCJS <- function(R,T,k,phi,p,control=NULL,hess=FALSE){
    ## Compute initial values on logistic scale
    inits <- pmax(pmin(c(logit$linkfun(phi),
                         logit$linkfun(p)),
                       logit$linkfun(.999)),
                  logit$linkfun(.001))

    ## Maximize likelihood and compute Hessian
    optim(inits,llhdCJS,llhdCJSgr,R=R,T=T,k=k,
          method="BFGS",debug=FALSE,
          control=c(list(fnscale=-1,trace=10),control),
          hessian=hess)
}

mleCJSnlm <- function(R,T,k,phi,p,print.level=1,hess=FALSE){
    ## Compute initial values on logistic scale
    inits <- pmax(pmin(c(logit$linkfun(phi),
                         logit$linkfun(p)),
                       logit$linkfun(.999)),
                  logit$linkfun(.001))

    ## Maximize likelihood and compute Hessian
    stats::nlm(llhdCJScombined,inits,R=R,T=T,k=k,llhdscale=-1,
               stepmax=3,print.level=print.level,hessian=hess,
               check.analyticals=TRUE)
}

mleCJSnlmtest <- function(R,T,k,phi,p,print.level=1,hess=FALSE){
    ## Compute initial values on logistic scale
    inits <- pmax(pmin(c(logit$linkfun(phi),
                         logit$linkfun(p)),
                       logit$linkfun(.999)),
                  logit$linkfun(.001))

    ## Maximize likelihood and compute Hessian
    stats::nlm(llhdCJScombinedtest,inits,R=R,T=T,k=k,llhdscale=-1,
               print.level=print.level,hessian=hess,
               check.analyticals=TRUE)
}

mleCJSrootSolve <- function(R,T,k,phi,p){
    ## Compute initial values on logistic scale
    inits <- pmax(pmin(c(logit$linkfun(phi),
                         logit$linkfun(p)),
                       logit$linkfun(.999)),
                  logit$linkfun(.001))

    rootSolve::multiroot(llhdCJSgr,inits,R=R,T=T,k=k,jacfunc=llhdCJShess,jactype="fullusr")
}
