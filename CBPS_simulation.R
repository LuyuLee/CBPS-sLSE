library(CBPS)

source('gama.R')
source('simulation_data.R')

alpha = matrix(c(0, -1, 0.5, -0.25, -0.1))
beta = matrix(c(20, 2.74, 1.37, 1.37, 1.37))

treat_size = 1
N = c(400,800)
n = N[2]
##q is the order of polynomial which is used to estimate E(Y|PS)
## q = (0, 2, 1)
q = 0
nitem = 10000
weight = NULL
set.seed(12345)
# psgen = c('logit', 'probit', 'tree')
psgen = 'tree'

for (q in 1:2)
{
for (n in N)
{
  
  # pre setting
  Beta_prob1 = vector(mode="list", length=nitem)
  Beta_over1 = vector(mode="list", length=nitem)
  Beta_exact1 = vector(mode="list", length=nitem)
  
  Beta_logit = vector(mode="list", length=nitem)
  
  Beta_logit2 = vector(mode="list", length=nitem)
  
  Beta_prob2 = vector(mode="list", length=nitem)
  Beta_over2 = vector(mode="list", length=nitem)
  Beta_exact2 = vector(mode="list", length=nitem)
  
  Beta_prob1_weight = vector(mode="list", length=nitem)
  Beta_over1_weight = vector(mode="list", length=nitem)
  Beta_exact1_weight = vector(mode="list", length=nitem)
  
  Beta_logit_weight = vector(mode="list", length=nitem)
  
  Beta_logit2_weight = vector(mode="list", length=nitem)
  
  Beta_prob2_weight = vector(mode="list", length=nitem)
  Beta_over2_weight = vector(mode="list", length=nitem)
  Beta_exact2_weight = vector(mode="list", length=nitem)
  
  for (i in 1 : nitem)
  {
 
    print(i)
    data = cbps_simulation(n = n, alpha = alpha, beta = beta,
                           delta = treat_size, psgen = psgen)
    Y = data$Y
    treat = data$treat
    # X is observed covariances, Z is latent logit variances
    Z = data[1 : (length(beta) - 1)]
    X = data[(length(beta)) : (length(data[1,]) - 2)]
    
    ### setting
    # both logit
    alpm_r = glm(formula = treat~X1+X2+X3+X4, data = data,
               family = binomial(link = 'probit'))
    alpha_r = alpm_r$coefficients
    pi_prob_right = as.numeric(predict(alpm_r, newdata = data, type = 'response'))
    
    fit1_over = CBPS(treat~X1+X2+X3+X4, data = data, ATT = 1)
    pi_over_right = fit1_over$fitted.values
    
    fit2_exact = CBPS(treat~X1+X2+X3+X4, data = data, ATT = 1, method = 'exact')
    pi_exact_right = fit2_exact$fitted.values
    
    logitm = glm(formula = treat~X1+X2+X3+X4, data = data, family = binomial(link = 'logit'))
    alpha_logit = logitm$coefficients
    pi_logit = logitm$fitted.values
    
    
    logitm = glm(formula = treat~Z1+Z2+Z3+Z4, data = data, family = binomial(link = 'logit'))
    alpha_logit_w = logitm$coefficients
    pi_logit_w = 1/(1+exp(-predict(logitm, newdata = data)))
    
    
    gama_prob_right = getgama_cbps(Y=Y, x=pi_prob_right, q=q)
    gama_over_right = getgama_cbps(Y=Y, x=pi_over_right, q=q)
    gama_exact_right = getgama_cbps(Y=Y, x=pi_exact_right, q=q)
    gama_logit = getgama_cbps(Y=Y, x=pi_logit, q=q)
    
    # misspecified ps
    alpm_w = glm(formula = treat~Z1+Z2+Z3+Z4, data = data,
                 family = binomial(link = 'probit'))
    alpha_w = alpm_w$coefficients
    pi_prob_wrong = as.numeric(predict(alpm_w, newdata = data, type = 'response'))
    
    fit1_over = CBPS(treat~Z1+Z2+Z3+Z4, data = data, ATT = 1)
    pi_over_wrong = fit1_over$fitted.values
    
    fit2_exact = CBPS(treat~Z1+Z2+Z3+Z4, data = data, ATT = 1, method = 'exact')
    pi_exact_wrong = fit2_exact$fitted.values
    
    # misspecified Ey
    gama_prob_wrong = getgama_cbps(Y=Y, x=pi_prob_wrong, q=q)
    gama_over_wrong = getgama_cbps(Y=Y, x=pi_over_wrong, q=q)
    gama_exact_wrong = getgama_cbps(Y=Y, x=pi_exact_wrong, q=q)
    gama_logit_wrong = getgama_cbps(Y=Y, x=pi_logit_w, q=q)
    
    # first scenario
    xd_prob = treat - pi_prob_right
    xd_cbps_over = treat - pi_over_right
    xd_cbps_exact = treat - pi_exact_right
    xd_logit = treat - pi_logit
    
    xa_prob_right = as.matrix(cbind(rep(1,n),Z)) %*% alpha_r
    
    # get beta hat value
    beta_prob1 = getpsr(Y = Y, dy = pi_prob_right, gama = gama_prob_right, X = xd_prob)
    beta_over1 = getpsr(Y = Y, dy = pi_over_right, gama = gama_over_right, X = xd_cbps_over)
    beta_exact1 = getpsr(Y = Y, dy = pi_exact_right, gama = gama_exact_right, X = xd_cbps_exact)
    beta_logit = getpsr(Y = Y, dy = pi_logit, gama = gama_logit, X = xd_logit)
    
    beta_prob1_weight = getwpsr(Y = Y, dy = pi_prob_right, gama = gama_prob_right, X = xd_prob, pi = pi_prob_right)
    beta_over1_weight = getwpsr(Y = Y, dy = pi_over_right, gama = gama_over_right, X = xd_cbps_over, pi = pi_over_right)
    beta_exact1_weight = getwpsr(Y = Y, dy = pi_exact_right, gama = gama_exact_right, X = xd_cbps_exact, pi = pi_exact_right)
    beta_logit_weight = getwpsr(Y = Y, dy = pi_logit, gama = gama_logit, X = xd_logit, pi = pi_logit)
    
    
    
    # 2nd scenario misspecified pi
    xd_prob = treat - pi_prob_wrong
    xd_cbps_over = treat - pi_over_wrong
    xd_cbps_exact = treat - pi_exact_wrong
    xd_logit_w = treat - pi_logit_w
    xa_prob_wrong = as.matrix(cbind(rep(1,n),X)) %*% alpha_w
    
    # get beta hat value
    beta_prob2 = getpsr(Y = Y, dy = pi_prob_wrong, gama = gama_prob_wrong, X = xd_prob)
    beta_over2 = getpsr(Y = Y, dy = pi_over_wrong, gama = gama_over_wrong, X = xd_cbps_over)
    beta_exact2 = getpsr(Y = Y, dy = pi_exact_wrong, gama = gama_exact_wrong, X = xd_cbps_exact)
    beta_logit2 = getpsr(Y = Y, dy = pi_logit_w, gama = gama_logit_wrong, X = xd_logit_w)
    
    beta_prob2_weight = getwpsr(Y = Y, dy = pi_prob_wrong, gama = gama_prob_wrong, X = xd_prob, pi = pi_prob_wrong)
    beta_over2_weight = getwpsr(Y = Y, dy = pi_over_wrong, gama = gama_over_wrong, X = xd_cbps_over, pi = pi_over_wrong)
    beta_exact2_weight = getwpsr(Y = Y, dy = pi_exact_wrong, gama = gama_exact_wrong, X = xd_cbps_exact, pi = pi_exact_wrong)
    beta_logit2_weight = getwpsr(Y = Y, dy = pi_logit_w, gama = gama_logit_wrong, X = xd_logit_w, pi = pi_logit_w)
    
    
    #save
    Beta_prob1[i] =  beta_prob1
    Beta_over1[i] = beta_over1
    Beta_exact1[i] = beta_exact1
    
    Beta_prob2[i] =  beta_prob2
    Beta_over2[i] = beta_over2
    Beta_exact2[i] = beta_exact2
    
    Beta_logit[i] = beta_logit
    Beta_logit2[i] = beta_logit2
    
    Beta_prob1_weight[i] =  beta_prob1_weight
    Beta_over1_weight[i] = beta_over1_weight
    Beta_exact1_weight[i] = beta_exact1_weight
    
    Beta_prob2_weight[i] =  beta_prob2_weight
    Beta_over2_weight[i] = beta_over2_weight
    Beta_exact2_weight[i] = beta_exact2_weight
    
    Beta_logit_weight[i] = beta_logit_weight
    Beta_logit2_weight[i] = beta_logit2_weight
  #save
  }  
  data_o = data.frame(Beta_logit = as.numeric(Beta_logit),
                      Beta_logit_weight = as.numeric(Beta_logit_weight),
                      
    
                      Beta_prob1 = as.numeric(Beta_prob1),
                      Beta_prob1_weight = as.numeric(Beta_prob1_weight),
                      
                      Beta_over1 = as.numeric(Beta_over1),
                      Beta_over1_weight = as.numeric(Beta_over1_weight),
                      
                      Beta_exact1 = as.numeric(Beta_exact1),
                      Beta_exact1_weight = as.numeric(Beta_exact1_weight),
                      
                      
                      Beta_logit2 = as.numeric(Beta_logit2),
                      Beta_logit2_weight = as.numeric(Beta_logit2_weight),
                      
                      Beta_prob2 = as.numeric(Beta_prob2),
                      Beta_prob2_weight = as.numeric(Beta_prob2_weight),
                      
                      Beta_over2 = as.numeric(Beta_over2),
                      Beta_over2_weight = as.numeric(Beta_over2_weight),
                      
                      
                      Beta_exact2 = as.numeric(Beta_exact2),
                      Beta_exact2_weight = as.numeric(Beta_exact2_weight)
                      
                      )
  
  means = data.frame(
                     Beta_logit = mean(as.numeric(Beta_logit)),
                     Beta_logit_weight = mean(as.numeric(Beta_logit_weight)),
                     Beta_prob1 = mean(as.numeric(Beta_prob1)),
                     Beta_prob1_weight= mean(as.numeric(Beta_prob1_weight)),
                     
                     Beta_over1 = mean(as.numeric(Beta_over1)),
                     Beta_over1_weight = mean(as.numeric(Beta_over1_weight)),
                     Beta_exact1 = mean(as.numeric(Beta_exact1)),
                     Beta_exact1_weight = mean(as.numeric(Beta_exact1_weight)),
                     
                     Beta_logit2 = mean(as.numeric(Beta_logit2)),
                     Beta_logit2_weight = mean(as.numeric(Beta_logit2_weight)),
                     Beta_prob2 = mean(as.numeric(Beta_prob2)),
                     Beta_prob2_weight = mean(as.numeric(Beta_prob2_weight)),
                     Beta_over2 = mean(as.numeric(Beta_over2)),
                     Beta_over2_weight = mean(as.numeric(Beta_over2_weight)),
                     Beta_exact2 = mean(as.numeric(Beta_exact2)),
                     Beta_exact2_weight = mean(as.numeric(Beta_exact2_weight)))
  
  data1 = data_o
  data2 = data_o - treat_size
  for (j in 1 : length(data1[1,]))
  {
    data1[,j] = data1[,j] - as.numeric(means[j])  
  }
  standdev = data.frame(Stand_DEV = diag((t(as.matrix(data1)) %*% as.matrix(data1)) ** 0.5) / (nitem ** 0.5))
  Rmse = data.frame(RMSE = diag((t(as.matrix(data2)) %*% as.matrix(data2)) ** 0.5) / (nitem ** 0.5) )
  ans = data.frame(abs_bias = t(as.matrix(abs(means - treat_size) )), standdev, Rmse)
  filename = paste(psgen,'_CBPS_', as.character(nitem), '_gama_q=', as.character(q),
                   '_num=', as.character(n),'.csv',sep = '')
  View(ans)
  write.csv(ans, filename)
}
}
