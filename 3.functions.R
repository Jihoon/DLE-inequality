#'Derive reference Pareto pdf from existing statistics (reference != base year distr)
#'@param xval vector of x variable (domain)
#'@param gini existing gini from data for the base year
#'@param avg  existing mean from data for the base year
#'@param min  lowest value of the existing distribution (if known)

RefPareto <- function(xval, gini, ave, min=0) {
  k = seq(0.01, 2, 0.01)
  
  base_gini = gini * ave / (ave - (1-k)*min)
  
  alpha = 0.5*(1/base_gini + 1)
  
  base_mean = (avg - (1-k)*min) / k
  base_mean2 = alpha * min / (alpha-1)
  
  sol = which((base_mean > base_mean2) ==FALSE)[1]
  
  k_sol = k[sol]
  alpha_sol = alpha[sol]
  d_sol = (1-k_sol) * min
  
  x0 = dpareto(xval, min, shape=alpha_sol)   # Gini = 1/(2*alpha-1)= 0.54 Staus quo
  
  return(x0)
}


#'Derive reference lognormal pdf from existing statistics 
#'@param xval vector of x variable (domain)
#'@param gini existing gini from data for the base year
#'@param avg  existing mean from data for the base year
#'@param min  lowest value of the existing distribution (if known)

RefLognorm <- function(xval, gini, avg, min=0){
  gini_baselognormal = gini *avg / (avg-min) # Assume scaler=1 for this base transformation
  
  # Parameters of the base log normal pdf
  sdlog_hat = qnorm((gini_baselognormal+1)/2)*sqrt(2) 
  mlog_hat = log(avg-min) - (sdlog_hat^2)/2 
  
  x0 = dlnorm(xval, mlog_hat, sdlog_hat)
  # F0 = plnorm(xval, mlog_hat, sdlog_hat)
  
  return(x0)
  # return(list(x0, F0))
}

DrawRefLognorm <- function(n, gini, avg, min=0){
  gini_baselognormal = gini *avg / (avg-min) # Assume scaler=1 for this base transformation
  
  # Parameters of the base log normal pdf
  sdlog_hat = qnorm((gini_baselognormal+1)/2)*sqrt(2) 
  mlog_hat = log(avg_ene-actual.min) - (sdlog_hat^2)/2 
  
  x = rlnorm(n, mlog_hat, sdlog_hat)
  
  return(x)
}

#'Transform distribution to achieve a target Gini
#'This translates the equations in the manuscript.
#'@param gini.base Gini of base year distribution
#'@param gini.tgt  Target Gini 
#'@param avg.base  Average of base year distribution
#'@param min.base  Minimum of base year distribution
#'@param min.tgt   Target minimum (e.g. DLE)
#'@param n_yr      Number of years until achievement

TransformDistr <- function(gini.base, gini.tgt, avg.base, min.base=0, min.tgt, n_yr) {
  gini.r = gini.tgt/gini.base
  
  scaler.k = gini.r*min.tgt/(avg.base*(1-gini.r) + gini.r*min.base)  # k derived from Eqn (8)
  
  d = min.tgt - scaler.k*min.base # D = min.tgt in Eqn (7)
  # d = D - scaler.k*min.base
  
  avg.new = avg.base*scaler.k + d
  
  growth.yr = (avg.new/avg.base)^(1/n_yr)-1 # Overall growth rate per year per capita (r-p), Eqn (9)
  
  return(list(k=scaler.k, d=d, avg.new=avg.new, gr=growth.yr))
}


TopShare2FillBottom <- function(tot.bottom, distr, gini.base, avg.base, min.base) {
  k = seq(0.01, 2, 0.01)
  
  base_gini = gini * ave / (ave - (1-k)*min)
  
  alpha = 0.5*(1/base_gini + 1)
  
  base_mean = (avg - (1-k)*min) / k
  base_mean2 = alpha * min / (alpha-1)
  
  sol = which((base_mean > base_mean2) ==FALSE)[1]
  
  k_sol = k[sol]
  alpha_sol = alpha[sol]
  d_sol = (1-k_sol) * min
  
  x0 = dpareto(xval, min, shape=alpha_sol)   # Gini = 1/(2*alpha-1)= 0.54 Staus quo
  
  return(x0)
}
