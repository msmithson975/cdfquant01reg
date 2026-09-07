##### Arcsinh_Arcsinh_outer_W_3 #####
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

Vtransf <- function(x, theta){
  return ((exp(theta)*x)/(1-x+exp(theta)*x))}

Vinv <- function(x, theta){
  return (x/(exp(theta) + x - x*exp(theta)))}

Wtransf <- function(x, theta){
  return (sinh(asinh(x) + theta))}

Winv <- function(x, theta){
  return (sinh(asinh(x) - theta))}

Utransf <- function(x, mu, sigma){
  return ((x-mu)/sigma)}

Uinv <- function(x, mu, sigma){
  return (mu+x*sigma)}

CDF <- function(x, mu, sigma, theta){
  if(x==0){
    return (sigma*exp(theta))}
  if(x==1){
    return (sigma*exp(-theta))}
  else{
    return (Arcsinhf(Wtransf(Utransf(InvArcsinhf(x),mu,sigma),theta)))}
}

# density function
dnsity <- function(x, mu, sigma, theta){
  if(x==0){
    return (sigma*exp(theta))}
  if(x==1){
    return (sigma*exp(-mu))}
  else{
    return ((CDF(x + 0.000001, mu, sigma, theta)-CDF(x, mu, sigma, theta))/0.000001)}
}

inverse_CDF <- function(x, mu, sigma, theta){
  return(Arcsinhf(Uinv(Winv(InvArcsinhf(x),theta),mu,sigma)))
}

# random draws function
rsamp <- function(n, mu, sigma, theta){
  return (inverse_CDF(runif(n, min = 0, max = 1),mu,sigma,theta)) }

# log_lik function
log_lik_Arcsinh_Arcsinh_outer_W_3 <- function(i, prep) {
  mu <- brms::get_dpar(prep, "mu", i = i)
  sigma <- brms::get_dpar(prep, "sigma", i = i)
  theta <- brms::get_dpar(prep, "theta", i = i)
  y <- prep$data$Y[i]
  log(dnsity(y, mu, sigma, theta))
}

# posterior_predict function
posterior_predict_Arcsinh_Arcsinh_outer_W_3 <- function(i, prep, ...) {
  mu <- brms::get_dpar(prep, "mu", i = i)
  sigma <- brms::get_dpar(prep, "sigma", i = i)
  theta <- brms::get_dpar(prep, "theta", i = i)
  ndraws <- prep$ndraws
  rsamp(ndraws, mu, sigma, theta)
}

# posterior_epred function
posterior_epred_Arcsinh_Arcsinh_outer_W_3 <- function(prep) {
  mu <- brms::get_dpar(prep, "mu")
  sigma <- brms::get_dpar(prep, "sigma")
  theta <- brms::get_dpar(prep, "theta")
  inverse_CDF(0.5,mu,sigma,theta)
}

# This defines the custom family.
Arcsinh_Arcsinh_outer_W_3 <- custom_family(
  "Arcsinh_Arcsinh_outer_W_3", dpars = c("mu", "sigma", "theta"),
  links = c("identity", "log", "identity"),
  lb = c(NA, 0, NA), ub = c(NA, NA, NA),
  type = "real"
)
# This is the piecewise log-likelihood. 
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
  
  real InvArcsinh(real x){
    return (1-2*x)/(2*x*(x-1));}
    
  real WTransform(real x, real theta){
    return sinh(asinh(x) + theta);}
    
  real U(real x, real mu, real sigma){
    return (x-mu)/sigma;}
  
  real CDFF(real x, real mu, real sigma, real theta){
  if(x==0){
    return sigma*exp(theta);}
  if(x==1){
    return sigma*exp(-theta);}
  else{
    return Arcsinh(WTransform(U(InvArcsinh(x),mu,sigma),theta));}
}
  
  real Arcsinh_Arcsinh_outer_W_3_lpdf(real x, real mu, real sigma, real theta){
    if(x==0){
      return log(sigma*exp(theta));}
    if(x==1){
      return log(sigma*exp(-theta));}
    else{
       return log((CDFF(x + 0.000001, mu, sigma, theta)-CDFF(x, mu, sigma, theta))/0.000001);}
}'
# We make an object that tells brms that the object stan_funs contains 
# STAN-code from the functions block of the STAN program. 
my_family <- Arcsinh_Arcsinh_outer_W_3
stanvars <- stanvar(scode = stan_funs, block = "functions")
#
