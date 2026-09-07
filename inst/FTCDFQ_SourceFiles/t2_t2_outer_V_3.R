##### t2_t2_outer_V_3 #####
# Helper functions

T2f <- function(x) {
  return ((x/(2*sqrt((x)^2+2))) + 0.5)}

signumf <- function(x){
  if(x < 0){
    return (-1)}
  else{ 
    return (1)}
}

InvT2f <- function(x) {
  return (signumf(x-0.5)*(sqrt((1-2*x)^2))/(sqrt(2)*sqrt((1-x)*x)))}

InvT2vec <- function(x) {
  return (ifelse(x < 0.5, -(sqrt((1-2*x)^2))/(sqrt(2)*sqrt((1-x)*x)), (sqrt((1-2*x)^2))/(sqrt(2)*sqrt((1-x)*x))))
}

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
  if (x==0){
    return((sigma^2)*exp(theta)) 
  }
  if (x==1){
    return((sigma^2)*exp(-theta)) 
  }
  else{
    return(Vtransf(T2f(Utransf(InvT2vec(x),mu,sigma)),theta))
  }
}

# density function
dnsity <- function(x, mu, sigma, theta){
  if(x==0){
    return ((sigma^2)*exp(theta))}
  if (x==1){
    return((sigma^2)*exp(-theta))}
  else{
    return ((CDF(x + 0.000001, mu, sigma, theta)-CDF(x, mu, sigma, theta))/0.000001)}
}

inverse_CDF <- function(x, mu, sigma, theta){
  return(T2f(Uinv(InvT2vec(Vinv(x,theta)),mu,sigma)))
}

# random draws function
rsamp <- function(n, mu, sigma, theta){
  return(T2f(Uinv(InvT2vec(Vinv(runif(n, min = 0, max = 1),theta)),mu,sigma)))}

# log_lik function
log_lik_t2_t2_outer_V_3 <- function(i, prep) {
  mu <- brms::get_dpar(prep, "mu", i = i)
  sigma <- brms::get_dpar(prep, "sigma", i = i)
  theta <- brms::get_dpar(prep, "theta", i = i)
  y <- prep$data$Y[i]
  log(dnsity(y, mu, sigma, theta))
}

# posterior_predict function
posterior_predict_t2_t2_outer_V_3 <- function(i, prep, ...) {
  mu <- brms::get_dpar(prep, "mu", i = i)
  sigma <- brms::get_dpar(prep, "sigma", i = i)
  theta <- brms::get_dpar(prep, "theta", i = i)
  ndraws <- prep$ndraws
  rsamp(ndraws, mu, sigma, theta)
}

# posterior_epred function
posterior_epred_t2_t2_outer_V_3 <- function(prep) {
  mu <- brms::get_dpar(prep, "mu")
  sigma <- brms::get_dpar(prep, "sigma")
  theta <- brms::get_dpar(prep, "theta")
  inverse_CDF(0.5,mu,sigma,theta)
}

# This defines the custom family.
t2_t2_outer_V_3 <- custom_family(
  "t2_t2_outer_V_3", dpars = c("mu", "sigma", "theta"),
  links = c("identity", "log", "identity"),
  lb = c(NA, 0, NA), ub = c(NA, NA, NA),
  type = "real"
)
# This is the piecewise log-likelihood. 
stan_funs <- '
  real T2(real x){
    return (x/(2*sqrt((x)^2+2))) + 0.5;}
    
  real signum(real x) {
    if(x < 0){
      return -1;}
    if(x > 0){
      return 1;}
    else{ 
      return 0;}
    }
  
  real InvT2(real x){
    return signum(x-0.5)*(sqrt((1-2*x)^2))/(sqrt(2)*sqrt((1-x)*x));}
    
  real VTransform(real x, real theta){
    return (exp(theta)*x)/(1-x+exp(theta)*x);}
    
  real U(real x, real mu, real sigma){
    return (x-mu)/sigma;}
  
    real CDFF(real x, real mu, real sigma, real theta){
    if (x==0){return log(pow(sigma,2)*exp(theta));
    }
    if (x==1){return log(pow(sigma,2)*exp(-theta));
    }
    else{ 
    return VTransform(T2(U(InvT2(x),mu,sigma)),theta);
    }
  }
  
  real t2_t2_outer_V_3_lpdf(real x, real mu, real sigma, real theta){
    if(x==0){
      return log(pow(sigma,2)*exp(theta));}
    if(x==1){
      return log(pow(sigma,2)*exp(-theta));}
    else{
       return log((CDFF(x + 0.000001, mu, sigma, theta)-CDFF(x, mu, sigma, theta))/0.000001);}
}'
# We make an object that tells brms that the object stan_funs contains 
# STAN-code from the functions block of the STAN program. 
my_family <- t2_t2_outer_V_3
stanvars <- stanvar(scode = stan_funs, block = "functions")
#
