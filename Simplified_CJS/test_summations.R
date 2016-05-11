## Load code
source("../TruncatedCJS/R/utilities.R")

## Recovery probability function
compP = function(phi,p,k){
  phi^(1:k)*(1-p)^(0:(k-1))*p
}

## Parameters
phi = .8
p = .2
k = 10

## Compute cell probabilities
P = compP(phi,p,k)

## 1) \sum_{j=1}^k P_j

sum(P) - 
  p*phi*(1-phi^k*(1-p)^k)/(1-phi*(1-p))

## 2) \sum_{j=1}^k jP_j/phi

sum((1:k)*P/phi) -
    p * sum((1:k) * ((1-p)*phi)^(1:k-1))

## Raw maxima output
sum((1:k)*P/phi) -
    ((p*(1-(1-p)^k*phi^k))/(1-(1-p)*phi)-((p-1)*p*phi*(1-(1-p)^k*phi^k))/(1-(1-p)*phi)^2-(k*(1-p)^k*p*phi^k)/(1-(1-p)*phi))

## Simplified
sum((1:k)*P/phi) -
    (p*(1 -(k+1)*phi^k*(1-p)^k+k*phi^(k+1)*(1-p)^(k+1))/(1-phi*(1-p))^2)

## 3) \sum_{j=2}^k j(j+1)p_j/\phi^2

## Raw maxima output
sum((1:k)*(1:k-1)*P/phi^2) -
    (-(2*(p-1)*p*(1-(1-p)^k*phi^k))/(1-(1-p)*phi)^2+(2*(p-1)^2*p*phi*(1-(1-p)^k*phi^k))/(1-(1-p)*phi)^3-(k^2*(1-p)^k*p*phi^(k-1))/(1-(1-p)*phi)-(k*(1-p)^k*p*phi^(k-1))/(1-(1-p)*phi)+(2*k*(1-p)^k*(p-1)*p*phi^k)/(1-(1-p)*phi)^2)



