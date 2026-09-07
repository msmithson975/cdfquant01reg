##### t2_t2_ inner_V #####
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

# density function
dnsity <- function(x, mu, sigma, u){
  dex <- ifelse(x==0, T2f(Utransf(InvT2vec(Vtransf(u/(2*u+1),mu)),sigma)), 
                ifelse(x==1, 1 - T2f(Utransf(InvT2vec(Vtransf((u+1)/(2*u+1),mu)),sigma)), 
                       (T2f(Utransf(InvT2vec(Vtransf((u+x)/(2*u+1) + 0.000001,mu)),sigma)) - T2f(Utransf(InvT2vec(Vtransf((u+x)/(2*u+1),mu)),sigma)))/(0.000001*(2*u + 1))))
  return(dex)
}

CDF <- function(x, mu, sigma, u){
  ifelse(x < -u,0, ifelse(x > 1+u,1, T2f(Utransf(InvT2vec(Vtransf((u+x)/(2*u+1),mu)),sigma))))}

inverse_CDF <- function(x, mu, sigma, u){
  return((2*u+1)*(Vinv(T2f(Uinv(InvT2f(x),sigma)),mu))-u)
}

# random draws function
rsamp <- function(n, mu, sigma, u){
  return((2*u+1)*(Vinv(T2f(Uinv(InvT2vec(runif(n, min = 0, max = 1)),sigma)),mu))-u) }

# log_lik function
log_lik_Ex_t2_t2_inner_V <- function(i, prep) {
  mu <- brms::get_dpar(prep, "mu", i = i)
  sigma <- brms::get_dpar(prep, "sigma", i = i)
  u <- brms::get_dpar(prep, "u", i = i)
  y <- prep$data$Y[i]
  log(dnsity(y, mu, sigma, u))
}

# posterior_predict function
posterior_predict_Ex_t2_t2_inner_V <- function(i, prep, ...) {
  mu <- brms::get_dpar(prep, "mu", i = i)
  sigma <- brms::get_dpar(prep, "sigma", i = i)
  u <- brms::get_dpar(prep, "u", i = i)
  ndraws <- prep$ndraws
  rsamp(ndraws, mu, sigma, u)
}

# posterior_epred function
posterior_epred_Ex_t2_t2_inner_V <- function(prep) {
  mu <- brms::get_dpar(prep, "mu")
  sigma <- brms::get_dpar(prep, "sigma")
  u <- brms::get_dpar(prep, "u")
  inverse_CDF(0.5,mu,sigma,u)
}

# This defines the custom family.
Ex_t2_t2_inner_V <- custom_family(
  "Ex_t2_t2_inner_V", dpars = c("mu", "sigma", "u"),
  links = c("identity", "log", "log"),
  lb = c(NA, 0, 0), ub = c(NA, NA, NA),
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
  
  real Ex_t2_t2_inner_V_lpdf(real x, real mu, real sigma, real u){
    if(x==0){
      return log(T2(U(InvT2(VTransform(u/(2*u+1),mu)),sigma)));
    }
    if(x==1){
      return log(1 - T2(U(InvT2(VTransform((u+1)/(2*u+1),mu)),sigma)));
    }
    else{
       return log(T2(U(InvT2(VTransform((u+x)/(2*u+1)+0.00001,mu)),sigma)) - T2(U(InvT2(VTransform((u+x)/(2*u+1),mu)),sigma))) - 
    log(2*u + 1);}
}'
# We make an object that tells brms that the object stan_funs contains 
# STAN-code from the functions block of the STAN program. 
my_family <- Ex_t2_t2_inner_V
stanvars <- stanvar(scode = stan_funs, block = "functions")
#
