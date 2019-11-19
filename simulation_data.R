simulation_data = function(n = 500, rownum = 6,
                           gamalpha = 1, alpha = alpha, beta = beta,
                           delta = 0.75)
{
  #### X generate   1 - 6  X1X2X3 continuous  456 binary
  library(MASS)
  # N = c(500, 2000)
  # n = N[1]
  sigma = rep(0.5, rownum * rownum)
  covariance <- matrix(sigma,nrow=rownum,byrow=TRUE)
  for (i in 1 : rownum)
    covariance[i, i] = 1
  data <- mvrnorm(n=n, mu=rep(0, rownum), covariance)
  X = data.frame(data)
  X[, (rownum/2 + 1):rownum] = ifelse(X[, (rownum/2 + 1):rownum]<0, 1, 0)
  x = as.matrix(cbind(rep(1, n), X))
  
  #############generate PS by logit model
  ps = 1 / (exp(-(x %*% alpha)) + 1)
  D = rbinom(n, 1, ps)
  summary(D)
  
  #############generate Y
  
  Y = x %*% beta + D * delta + rnorm(n)
  summary(Y)
  return(data.frame(X0 = 1, X, Y = Y, treat = D))
}


cbps_simulation = function(n = 500, 
                           alpha = alpha, beta = beta,
                           delta = 1, psgen = 'logit', trsd = c(-0.8, 0.3, 0, 0, 0))
{
  # psgen = c('logit', 'probit', 'tree')
  # Z are observed covariances and non-linear transformed from X which are variance of true model
  # set.seed(11545)
  library(MASS)
  covariance = diag(length(beta) - 1)  
  data <- mvrnorm(n=n, mu=rep(0, length(beta) - 1), covariance)
  X = data.frame(data)
  x = as.matrix(cbind(rep(1, n), X))
  attach(X)
  Z1 = exp(X1/2)
  Z2 = X2/(1 + exp(X1)) + 10
  Z3 = (X1 * X3 / 25 + 0.6)**3
  Z4 = (X1 + X4 + 20)**2
  detach(X)
  Z = data.frame(Z1, Z2, Z3, Z4)
  
  ## get Y_true and logit(pi)
  if (psgen == 'logit'){
    ps = 1 / (exp(-(x %*% alpha)) + 1)
    D = rbinom(n, 1, ps)
  }
  if (psgen == 'probit'){
    ps = pnorm(x %*% alpha)
    D = rbinom(n, 1, ps)
  }
  if (psgen == 'tree')
    D = tree_ps(X, trsd = trsd, bnprob = 0.8)
  summary(D)
  
  Y = x %*% beta + D * delta + rnorm(n)
  summary(Y)
  return(data.frame(Z, X, Y = Y, treat = D))
}

cbps_censor_simulation = function(n = 500, censorrate = 0.5,
                           alpha = alpha, beta = beta,
                           delta = 210)
{
  # Z are observed covariances and non-linear transformed from X which are variance of true model
  library(MASS)
  source('survival_data.R')
  covariance = diag(length(beta) - 1)  
  data <- mvrnorm(n=n, mu=rep(0, length(beta) - 1), covariance)
  X = data.frame(data)
  x = as.matrix(cbind(rep(1, n), X))
  attach(X)
  Z1 = exp(X1/2)
  Z2 = X2/(1 + exp(X1)) + 10
  Z3 = (X1 * X3 / 25 + 0.6)**3
  Z4 = (X1 + X4 + 20)**2
  detach(X)
  Z = data.frame(Z1, Z2, Z3, Z4)
  
  ## get Y_true and logit(pi)
  ps = 1 / (exp(-(x %*% alpha)) + 1)
  D = rbinom(n, 1, ps)
  summary(D)
  
  Y = x %*% beta + D * delta + rnorm(n)
  summary(Y)
  censordf = survival_data(y = Y, censorrate = censorrate)
  c = censordf$c
  Delta = censordf$Delta
  return(data.frame(Z, X, Y = Y, treat = D, y_cens = c, Delta))
}


# trsd = c(-0.8, 0.3, 0, 0, 0)

tree_ps = function(X, trsd = rep(-1.058, 5), bnprob = 0.8)
{
  # trsd is thredshold  t1 t2 t3 t4 t5
  t1 = trsd[1]
  t2 = trsd[2]
  t3 = trsd[3]
  t4 = trsd[4]
  t5 = trsd[5]
  set.seed(12345)
  U = rnorm(length(X[,1]), sd = 0.15)
  ##rule##
  treat = rep(0, length(X[,1]))
  tot0 = 0
  tot1 = 0 
  tot2 = 0 
  tot3 = 0
  attach(X)
  for (i in 1 : length(X[,1]))
  {
    if ((X1[i] < t1) | (X2[i] < t1) | (X3[i] < t1))
    {
      if (X4[i] < t2)
      {
        tot0 = tot0 + 1
        next
      }
    }
    if (X4[i] > t3)
    {
      treat[i] = rbinom(1, 1, prob = (bnprob + U[i]))
      if ((bnprob + U[i]) > 1)
        treat[i] = 1
      tot1 = tot1 + 1
      next
    }
    if (X2[i] > t4)
    {
      treat[i] = rbinom(1, 1, prob = bnprob + U[i])
      if ((bnprob + U[i]) > 1)
        treat[i] = 1
      tot2 = tot2 + 1
      next
    }
    if (X3[i] > t5)
    {
      treat[i] = rbinom(1, 1, prob = bnprob + U[i])
      if ((bnprob + U[i]) > 1)
        treat[i] = 1
      tot3 = tot3 + 1
      next
    }
  }
  detach(X)
  summary(treat)
  return(treat)
}