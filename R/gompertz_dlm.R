#### FUNCAO SIR
### FALTOU OS Y RESTANTES DOS DADOS

sir_gompertz_dlm <- function(fit, data, resampling_size = nrow(fit$mu)){

  ## Gompertz: gomp = a*exp(b*x)

    y = data$y
    sigma = sqrt(median(fit$sig2))

    likelihood <- function(par){
      gomp = par[1]*exp(par[2]*data$x)     #### dlm usa normalidade
      prod(dnorm(y, mean = log(gomp), sd = sigma))
    }

  ### Método SIR

  ## sampling A e B
  A = rbeta(500, 1, 10000)
  B = rbeta(500, 1, 10)

  ## Distribuição conjunta
  df = data.frame(A = rep(A, each = 500), B = rep(B, times = 500))

  ## Assumindo prioris uniformes, logo, distribuição a posteriori é proporcional a verossimilhança
  post.dist <- apply(df, 1, likelihood)

  pesos = post.dist/(dbeta(df[,1], 1, 10000)*dbeta(df[,2], 1, 10))
  probs = pesos/sum(pesos)

  res = sample(500*500, resampling_size, replace = T, prob = probs)
  res_A = df[res,1]
  res_B = df[res,2]

  return(data.frame(A = res_A, B = res_B))

}
