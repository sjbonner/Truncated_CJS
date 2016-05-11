## Packages
library(numDeriv)
library(ggplot2)
library(dplyr)

## Source code
source("TruncatedCJS/R/llhd_cjs.R")
source("TruncatedCJS/R/utilities.R")

## Parameters
T = 10 # Occasions

p = rep(.9, T - 1) # Capture probabilities
phi = rep(.1, T - 1) # Survival probabilities

eta = c(log(phi) - log(1 - phi), log(p) - log(1 - p)) ## Parameters on logit scale

n = rep(10000, T - 1) # Releases

V = sapply(2:(T - 1), function(k) {
  cat(k, "...\n")
  
  ## Compute cell probabilities
  P = PCJS(phi, p, T, k)
  
  ## Compute expected data
  EM = n * P
  
  ## Compute Hessian matrix numerically
  R = cbind(n, EM)
  H = hessian(llhdCJS,
              eta,
              R = R,
              T = T,
              k = k)
  
  ## Compute Hessian matrix analytically (Not working)
  ## H = llhdCJShess(c(phi,p),R,T,k)
  
  ## Extract variance estimates
  var = diag(-solve(H))
})

# Vdf=data.frame(Variance=as.vector(V),
#                Parameter=factor(rep(rep(c("phi","p"),c(T-1,T-1)),T-2)),
#                Index=factor(rep(1:(T-1),2*(T-2))),
#                Occasion=factor(rep(c(1:(T-1),2:T),T-2)),
#                k=rep(2:(T-1),rep(2*(T-1),T-2)))
#
# print(ggplot(subset(Vdf,Index != T-1),aes(x=k,y=Variance,group=Occasion)) +
#         geom_line() + geom_point(aes(color=Occasion)) + facet_grid(Parameter ~ .))

## Convert to a data frame with variance as proportion of full data value
Vprop = t(apply(V, 1, function(v)
  v / tail(v, 1)))

Vpropdf = data.frame(
  Variance = as.vector(Vprop),
  Parameter = factor(rep(rep(
    c("phi", "p"), c(T - 1, T - 1)
  ), T - 2)),
  Index = factor(rep(1:(T - 1), 2 * (T - 2))),
  Occasion = factor(rep(c(1:(
    T - 1
  ), 2:T), T - 2)),
  k = rep(2:(T - 1), rep(2 * (T - 1), T - 2))
)

print(
  ggplot(
    subset(Vpropdf, Index != T - 1),
    aes(x = k, y = Variance, group = Occasion)
  ) +
    geom_line() + 
    geom_point(aes(color = Occasion)) + 
    facet_grid(Parameter ~ .) +
    ggtitle(paste0("phi=", phi[1], ", p=", p[1], ", n=", n[1])) #+ coord_trans(y="log10")
)

## Compute maximum inflatioin by parameter and k
subset(Vpropdf,Index != T-1) %>% group_by(Parameter,k) %>% summarise(max(Variance))
