##### t2_t2_inner_V #####
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

# CDF function
CDF <- function(x, mu, sigma){
  dex <- ifelse(x==0, (sigma^2)*exp(mu), 
                ifelse(x==1, (sigma^2)*exp(-mu), 
                       (T2f(Utransf(InvT2vec(Vtransf(x,mu)),sigma)))))
  return(dex)
}

# density function
dnsity <- function(x, mu, sigma){
  dex <- ifelse(x==0, (sigma^2)*exp(mu), 
                ifelse(x==1, (sigma^2)*exp(-mu), 
                       (T2f(Utransf(InvT2vec(Vtransf(x + 0.000001,mu)),sigma)) - T2f(Utransf(InvT2vec(Vtransf(x,mu)),sigma)))/0.000001))
  return(dex)
}

inverse_CDF <- function(x, mu, sigma){
  return(Vinv(T2f(Uinv(InvT2f(x),sigma)),mu))
}

# random draws function
rsamp <- function(n, mu, sigma){
  return(Vinv(T2f(Uinv(InvT2vec(runif(n, min = 0, max = 1)),sigma)),mu)) }

# log_lik function
log_lik_t2_t2_inner_V <- function(i, prep) {
  mu <- brms::get_dpar(prep, "mu", i = i)
  sigma <- brms::get_dpar(prep, "sigma", i = i)
  y <- prep$data$Y[i]
  log(dnsity(y, mu, sigma))
}

# posterior_predict function
posterior_predict_t2_t2_inner_V <- function(i, prep, ...) {
  mu <- brms::get_dpar(prep, "mu", i = i)
  sigma <- brms::get_dpar(prep, "sigma", i = i)
  ndraws <- prep$ndraws
  rsamp(ndraws, mu, sigma)
}

# posterior_epred function
posterior_epred_t2_t2_inner_V <- function(prep) {
  mu <- brms::get_dpar(prep, "mu")
  sigma <- brms::get_dpar(prep, "sigma")
  inverse_CDF(0.5,mu,sigma)
}

# This defines the custom family.
t2_t2_inner_V <- custom_family(
  "t2_t2_inner_V", dpars = c("mu", "sigma"),
  links = c("identity", "log"),
  lb = c(NA, 0), ub = c(NA, NA),
  type = "real"
)
# This is the piecewise log-likelihood. 
# STAN lacks native cotangent or cosecant functions so these must be defined. 
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
    
  real VTransform(real x, real mu){
    return (exp(mu)*x)/(1-x+exp(mu)*x);}
    
  real U(real x, real sigma){
    return x/sigma;}
  
  real t2_t2_inner_V_lpdf(real x, real mu, real sigma){
    if(x==0){
      return log(pow(sigma,2)*exp(mu));
    }
    if(x==1){
      return log(pow(sigma,2)*exp(-mu));
    }
    else{
       return log(T2(U(InvT2(VTransform(x+0.00001,mu)),sigma)) - T2(U(InvT2(VTransform(x,mu)),sigma)));}
}'
# We make an object that tells brms that the object stan_funs contains 
# STAN-code from the functions block of the STAN program. 
my_family <- t2_t2_inner_V
stanvars <- stanvar(scode = stan_funs, block = "functions")
#
