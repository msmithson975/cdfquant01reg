##### Arcsinh_Cauchy_outer_W #####
# Helper functions
Cauchyf <- function(x){
  return ((atan(x)/pi) + 0.5) }

InvCauchyf <- function(x){
  return (tan(pi*x - pi/2)) }

signum <- function(x){
  if(x < 0){
    return (-1)}
  if(x > 0){
    return (1)}
  else{ 
    return (0)}
}

Arcsinhf <- function(x){
  return (1/(exp(-asinh(x))+1)) }

InvArcsinhf <- function(x){
  return ((1-2*x)/(2*x*(x-1)))}

Vtransf <- function(x, mu){
  return ((exp(mu)*x)/(1-x+exp(mu)*x))}

Vinv <- function(x, mu){
  return (x/(exp(mu) + x - x*exp(mu)))}

Wtransf <- function(x, mu){
  return (sinh(asinh(x) + mu))}

Winv <- function(x, mu){
  return (sinh(asinh(x) - mu))}

Utransf <- function(x, sigma){
  return (x/sigma)}

Uinv <- function(x, sigma){
  return (x*sigma)}

# density function
dnsity <- function(x, mu, sigma, u){
  if(x==0){
    return (Arcsinhf(Wtransf(Utransf(InvCauchyf(u/(2*u+1)),sigma),mu)))}
  if(x==1){
    return (1 - Arcsinhf(Wtransf(Utransf(InvCauchyf((u+1)/(2*u+1)),sigma),mu)))}
  else{
    return ((Arcsinhf(Wtransf(Utransf(InvCauchyf((u+x)/(2*u+1) + 0.000001),sigma),mu)) - Arcsinhf(Wtransf(Utransf(InvCauchyf((u+x)/(2*u+1)),sigma),mu)))/(0.000001*(2*u + 1)))}
}

CDF <- function(x, mu, sigma, u){
  ifelse(x < -u,0, ifelse(x > 1+u,1, Arcsinhf(Wtransf(Utransf(InvCauchyf((u+x)/(2*u+1)),sigma),mu))))}

inverse_CDF <- function(x, mu, sigma, u){
  return((2*u+1)*(Cauchyf(Uinv(Winv(InvArcsinhf(x),mu),sigma)))-u)
}

# random draws function
rsamp <- function(n, mu, sigma, u){
  return (inverse_CDF(runif(n, min = 0, max = 1),mu,sigma, u)) }

# log_lik function
log_lik_Ex_Arcsinh_Cauchy_outer_W <- function(i, prep) {
  mu <- brms::get_dpar(prep, "mu", i = i)
  sigma <- brms::get_dpar(prep, "sigma", i = i)
  u <- brms::get_dpar(prep, "u", i = i)
  y <- prep$data$Y[i]
  log(dnsity(y, mu, sigma, u))
}

# posterior_predict function
posterior_predict_Ex_Arcsinh_Cauchy_outer_W <- function(i, prep, ...) {
  mu <- brms::get_dpar(prep, "mu", i = i)
  sigma <- brms::get_dpar(prep, "sigma", i = i)
  u <- brms::get_dpar(prep, "u", i = i)
  ndraws <- prep$ndraws
  rsamp(ndraws, mu, sigma, u)
}

# posterior_epred function
posterior_epred_Ex_Arcsinh_Cauchy_outer_W <- function(prep) {
  mu <- brms::get_dpar(prep, "mu")
  sigma <- brms::get_dpar(prep, "sigma")
  u <- brms::get_dpar(prep, "u")
  inverse_CDF(0.5,mu,sigma,u)
}

# This defines the custom family.
Ex_Arcsinh_Cauchy_outer_W <- custom_family(
  "Ex_Arcsinh_Cauchy_outer_W", dpars = c("mu", "sigma", "u"),
  links = c("identity", "log", "log"),
  lb = c(NA, 0, 0), ub = c(NA, NA, NA),
  type = "real"
)
# This is the piecewise log-likelihood. 
# STAN lacks native cotangent or cosecant functions so these must be defined. 
stan_funs <- '
  real Arcsinh(real x){
    return 1 / (exp(-asinh(x))+1);}
    
  real signum(real x) {
    if(x < 0){
      return -1;}
    if(x > 0){
      return 1;}
    else{ 
      return 0;}
    }
  
  real InvCauchy(real x){
    return tan(pi()*x - pi()/2);}
    
  real WTransform(real x, real mu){
    return sinh(asinh(x) + mu);}
    
  real U(real x, real sigma){
    return x/sigma;}
  
  real Ex_Arcsinh_Cauchy_outer_W_lpdf(real x, real mu, real sigma, real u){
    if(x==0){
      return log(Arcsinh(WTransform(U(InvCauchy(u/(2*u+1)),sigma),mu)));
    }
    if(x==1){
      return log(1 - Arcsinh(WTransform(U(InvCauchy((u+1)/(2*u+1)),sigma),mu)));
    }
    else{
       return log(Arcsinh(WTransform(U(InvCauchy((u+x)/(2*u+1)+0.00001),sigma),mu)) - Arcsinh(WTransform(U(InvCauchy((u+x)/(2*u+1)),sigma),mu))) - 
    log(2*u + 1);}
}'
# We make an object that tells brms that the object stan_funs contains 
# STAN-code from the functions block of the STAN program. 
my_family <- Ex_Arcsinh_Cauchy_outer_W
stanvars <- stanvar(scode = stan_funs, block = "functions")
#
