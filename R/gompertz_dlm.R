#### SIR function

sir_gompertz_dlm <- function(fit, data, resampling_size = nrow(fit$mu)){

  ## Gompertz: gomp = a*exp(b*x)

  y = data$y  ## log(qx)
  sigma = median(sqrt(0.01*fit$sig2))

  likelihood <- function(par){
    gomp = par[1]*exp(par[2]*data$x)     #### normality in DLM
    prod(dnorm(y, mean = log(gomp), sd = sigma))
  }

  ### SIR method

  ## sampling A and B
  A = rbeta(500, 1, 10000)
  B = rbeta(500, 1, 10)

  ## Joint distribution
  df = data.frame(A = rep(A, each = 500), B = rep(B, times = 500))

  ## Assuming uniform priori, posteriori distribution is proportional to the likelihood
  post.dist <- apply(df, 1, likelihood)

  pesos = post.dist/(dbeta(df[,1], 1, 10000)*dbeta(df[,2], 1, 10))
  probs = pesos/sum(pesos)

  res = sample(500*500, resampling_size, replace = T, prob = probs)
  res_A = df[res,1]
  res_B = df[res,2]

  return(data.frame(A = res_A, B = res_B))

}
