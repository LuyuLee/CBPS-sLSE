getgama = function(Y=Y, X=X , alpha = alpha, q=0, weight = NULL)
{
  # creat var name
  if (q == 0)
  {
    return(mean(as.matrix(Y)))
  }
  f <- ""
  nextx <- ""
  x = as.matrix(cbind(list(rep(1, length(X[, 1]))), X)) %*% alpha
  dx = data.frame(Y)
  f <- ""
  nextx <- ""
  if (q>1) {
    for (ii in 1:(q-1)) {
      dx = data.frame(dx, x**ii)
      nextx <- paste("x",ii, sep="")
      if (ii==1) {name <- nextx}
      if (ii>1) {name <- c(name, nextx)}
      f <- paste(f, nextx, "+", sep="")
    }
    f <- paste(f, "x", ii+1, sep="")
  } else if (q==1) {
    f <- "x1"
  }
  dx = data.frame(dx, x**q)
  formulation = paste("Y","~", f, sep="")  
  name = ""
  for (ii in 1:q ) {
    nextx <- paste("x",ii, sep="")
    if (ii==1) {name <- nextx}
    if (ii>1) {name <- c(name, nextx)}
  }
  name = c("Y", name)
  names(dx) = name
  lmy = lm(formula = formulation, data = dx, weights = weight)
  gama = lmy$coefficients
  return(as.numeric(gama))
}


getgama_cbps = function(Y=Y, x=X, q=0, weight = NULL)
{
  # creat var name
  if (q == 0)
  {
    return(mean(as.matrix(Y)))
  }
  f <- ""
  nextx <- ""
  dx = data.frame(Y)
  f <- ""
  nextx <- ""
  if (q>1) {
    for (ii in 1:(q-1)) {
      dx = data.frame(dx, x**ii)
      nextx <- paste("x",ii, sep="")
      if (ii==1) {name <- nextx}
      if (ii>1) {name <- c(name, nextx)}
      f <- paste(f, nextx, "+", sep="")
    }
    f <- paste(f, "x", ii+1, sep="")
  } else if (q==1) {
    f <- "x1"
  }
  dx = data.frame(dx, x**q)
  formulation = paste("Y","~", f, sep="")  
  name = ""
  for (ii in 1:q ) {
    nextx <- paste("x",ii, sep="")
    if (ii==1) {name <- nextx}
    if (ii>1) {name <- c(name, nextx)}
  }
  name = c("Y", name)
  names(dx) = name
  lmy = lm(formula = formulation, data = dx, weights = weight)
  gama = lmy$coefficients
  return(as.numeric(gama))
}


getpsr = function(Y = Y, dy = Xa, gama = gama0, X = xd, weights = NULL)
{
  sum = gama[1]
  if (length(gama) > 1)
    for (i in 1 : (length(gama) - 1))
      sum = sum + gama[i + 1] * dy ** i
  y = as.matrix(Y - sum)
  lmy = lm(formula = y~X, weights = weights)
  beta = lmy$coefficients[2]
  return(beta)
}


getwpsr = function(Y = Y, dy = Xa, gama = gama0, X = xd, weights = NULL, pi = pi)
{
  if (is.null(weights))
  {
    df = data.frame(Y, dy, X, pi)
    df_wpsr = select_pi(df)
    Y = df_wpsr$Y
    pi = df_wpsr$pi
    dy = df_wpsr$dy
    X = df_wpsr$X
  }
  if (is.null(weights) == F)
  {
    df = data.frame(Y, dy, X, weights, pi)
    df_wpsr = select_pi(df)
    Y = df_wpsr$Y
    pi = df_wpsr$pi
    dy = df_wpsr$dy
    X = df_wpsr$X
    weights = df_wpsr$weights
  }
  IPW = ((1 - pi) * (pi)) ** 0.5
  sum = gama[1] 
  if (length(gama) > 1)
    for (i in 1 : (length(gama) - 1))
      sum = sum + gama[i + 1] * dy ** i
  y = as.matrix((Y - sum) / IPW)
  x = X / IPW
  lmy = lm(formula = y~x, weights = weights)
  beta = lmy$coefficients[2]
  return(beta)
} 
  
select_pi = function(df, upbound = 0.95, lowbound = 0.05)
{
  chose_num = list()
  item = 0
  for (i in 1: length(df$pi))
  {
    if (df$pi[i] > upbound)
      next
    else if (df$pi[i] < lowbound)
      next
    item = item + 1
    chose_num [item] = i
  }
  redf = df[as.numeric(chose_num),]
  return(redf)
}

