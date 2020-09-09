rm(list=ls())

M = 20000 # M is the number of hypothesis tests
L = 1000 # Number of tests
pi0 = 0.9 # proportion of null
beta_para = c(0.5,1)
uni_para = c(0,1)
thres = 0.1

# NOTE: ALL LIBRARIES ARE NOT NEEDED 
list.of.packages <- c("data.table", "fasttime",'plyr',"PerformanceAnalytics",
                      'imputeTS',"parallel","doParallel","doMC",'lubridate',
                      'anytime','xts','TTR','missRanger','quantmod')
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(list.of.packages, require, character.only = TRUE)


# LOAD ENVIRONMENT
if(Sys.info()['sysname'] == "Windows" ) {
  library(doParallel)
  registerDoParallel(cores=detectCores())
} else {
  library(doMC)
  registerDoMC(cores=detectCores())
}

gen_data = function(pi0,beta_para,uni_para) {
  pi1 = 1-pi0
  # y is the hidden variable
  # indicating null or non-null
  y = sample(0:1,M,replace = T,prob = c(pi0,pi1))
  mu = rep(0,M)
  mu[y==0] = runif(sum(y==0),min=uni_para[1],max=uni_para[2])
  mu[y==1] = rbeta(sum(y==1),shape1=beta_para[1],shape2=beta_para[2])
  mu
}
beta = function(p,para) {
  para[1]*(p^(para[1]-1)*(1-p)^(para[2]-1))
}
analytical = function(pi0,p,para) {
  pi1 = 1 - pi0
  (pi0) / (pi0 + pi1 * beta(p,para))
}
mu = gen_data(pi0,beta_para,uni_para)
post = analytical(pi0,mu,beta_para) 
post = 1 - post
hist(post)
count = get_pvalues(post,thres)


get_lfdr = function(pi0,beta_para,uni_para) {
  mu = gen_data(pi0,beta_para,uni_para)
  post = analytical(pi0,mu,beta_para) 
  post = 1 - post
  post
}

get_pvalues = function(post,thres) {
  ### Bonferroni correction
  bon = pmin(1,  length(post)* post)
  bon_thres = length(which(bon < thres))
  names(bon_thres) = 'bon'
  ### Benjamini-Hochberg Procedure
  lp <- length(post)
  i <- lp:1L
  o <- order(post, decreasing = TRUE)
  ro <- order(o)
  bh = pmin(1, cummin(lp/i * post[o]))[ro]
  bh_thres = length(which(bh < thres))
  names(bh_thres) = 'bh'
  c(bon_thres,bh_thres)
}

post = get_lfdr(pi0,beta_para,uni_para)

all_experiments = function(L,pi0,beta_para,uni_para) {
  grid.param <- expand.grid(1:L)
  fe <- foreach(param = iter(grid.param, by = "row"), 
                .verbose = TRUE, .errorhandling = "pass",  
                .multicombine = TRUE, .maxcombine = max(2, nrow(grid.param)),
                .export=c("analytical","beta","gen_data","get_lfdr","get_pvalues","M",'thres'),
    #names(Filter(is.function, mget(ls(".GlobalEnv")[-which(ls(".GlobalEnv") == 'all_experiments')])))
                .packages="foreach")
  
  fe$args <- fe$args[1]
  fe$argnames <- fe$argnames[1]
  
  results <- fe %dopar% {
    post = get_lfdr(pi0,beta_para,uni_para)
    count = get_pvalues(post,thres)
  }
  results = do.call(rbind,results)
  results
}

results = all_experiments(L,pi0,beta_para,uni_para)
print(apply(results,2,quantile))

stopImplicitCluster()


# Question 2

M = 20000
pi0 = 0.95


gen_data= function(M,pi0) {
  pi1 = 1-pi0 
  y = sample(0:1,M,replace = T,prob = c(pi0,pi1))
  mu = rep(0,M)
  mu[y==0] = rnorm(sum(y==0),0,0.5)
  mu[y==1] = rnorm(sum(y==1),2.5,0.5)
  mu
}

analytical = function(mu) {
  mu_ml = sum(mu) / length(mu)
  cov_ml = (1/length(mu)) * sum((mu - mu_ml) * (mu - mu_ml))
  1/(2 * pi *  abs(cov_ml))^0.5 * exp(-0.5 * (mu - mu_ml) * cov_ml^-1 * (mu - mu_ml))
}

convo = function(mu,z) {
  exp(-0.5 * (mu - z) * (var(mu) + var(z))^-1 * (mu - z)) / (2 * pi *  abs(var(mu) + var(z))^0.5)
}

conv_variance = function(pi0) {
  pi0* 0.5 + (1-pi0) * 0.5 + (pi0 * 0^2 + (1-pi0)*2.5^2 - (pi0*0+(1-pi0)*2.5)^2)
}

analytical_better = function(pi0,z) {
  pi0 * dnorm(z,0,sqrt(conv_variance(pi0))) + (1- pi0) * dnorm(z,2.5,sqrt(conv_variance(pi0)))
}


get_lfdr = function(pi0,mu) {
  pi1 = 1 - pi0
  (pi1 * dnorm(mu,2.5,sqrt(0.5))) / (pi0 * dnorm(mu,0,sqrt(0.5)) + pi1 * dnorm(mu,2.5,sqrt(0.5))) 
}

# E(mu | z)
posterior = function(pi0,z) {
  pi1 = 1 - pi0
  (pi1 * dnorm(z,2.5,0.5)) / (pi0 * dnorm(z,0,0.5) + pi1 * dnorm(z,2.5,0.5))
}

# Part A
## sanity check: should be equal to 0
mu = gen_data(M,pi0) 
print(var(mu) - conv_variance(pi0))
# Part B
z =  1 - get_lfdr(pi0,mu)
hist(z)
# Part C
post = 1 - posterior(pi0,mu)


# Question 3
### EM
M = 10000
alpha = 3^-1
iteration =100
for(i in 1:iteration) {
  mu = rnorm(M,0,sqrt(alpha))
  z =rnorm(M,mu,1)
  alpha_new = M / sum(1/(1+alpha) + z^2/(1+alpha)^2)
  alpha = alpha_new
}
print(sqrt(alpha))

## James Stein
alpha = 3^-1

js = function(draw) {
  l2 <- sum(draw^2)
  return((1 - (length(draw) - 2) / l2) *draw)
}
mu = rnorm(M,0,sqrt(alpha))
mu_js = js(mu)
var(mu_js)

# Hence, we conclude James Stein estimator have a larger variance

# Question 4 
### I did not have enough time I am affraid, I took 7courses this semester. My apologies.