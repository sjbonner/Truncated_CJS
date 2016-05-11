## Packages
library(maxLik)

source("../TruncatedCJS/R/utilities.R")

## Recovery probability function
compP = function(phi,p,k){
  phi^(1:k)*(1-p)^(0:(k-1))*p
}

## Likelihood function
llhd = function(phi,myp,mym){
  
 # browser()
  
  ##phi = logit$linkinv(eta)
  
  k = length(mym)-1
  
  P = compP(phi,myp,k)
  
  sum(mym[-(k+1)] * log(P)) + mym[(k+1)] * log(1- sum(P))
}

myvar = function(phi,p,n,k){
  
  P = compP(phi,p,k)
  
  phi^2/n/(sum((1:k)^2*P) + sum((1:k)*P)^2/(1-sum(P)))
  
}

## Simulate data
n = 10000
phi = .8
p = .2
k = 100


P = compP(phi,p,k)

m = rmultinom(1000,n,c(P,1-sum(P)))

## Maximize likelihood
ml = apply(m,2,function(m)
  maxLik(llhd,start=.5,myp=p,mym=m)$estimate
)

print(c(sd(ml),sqrt(myvar(phi,p,n,k))))

## Limit?
print(sqrt(-(-n*p*(1-p)/(phi*(1-phi*(1-p))^2) - n/((1-phi)*(1-(1-p)*phi))+(n*(p-1)^2*(1-phi))/(1-(1-p)*phi)^3)^-1))

print(sqrt((n*p/(1-phi*(1-p))^3 * (p/(1-phi) + (1 + phi*(1-p))/phi))^-1))
