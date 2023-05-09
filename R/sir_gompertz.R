
sir_gompertz <- function(fit, data, resampling_size = nrow(fit$post.samples$mcmc_theta)){

  ## Gompertz: gomp = a*exp(b*x)

  ## Verossimilhança
  if(fit$info$model == "binomial"){

    sim_Ex = rep(1000, length(data$Ex))
    sim_Dx = round(data$qx*sim_Ex)

    likelihood <- function(par){
      gomp = par[1]*exp(par[2]*data$x)
      q_x = 1 - exp(-gomp)
      prod(dbinom(sim_Dx, size = sim_Ex, prob = q_x))
    }

  }else if(fit$info$model == "poisson"){

    sim_Ex = rep(10000, length(data$Ex))
    sim_Dx = round(data$qx*sim_Ex)

    likelihood <- function(par){
      mx = par[1]*exp(par[2]*data$x)
      q_x = 1 - exp(-mx)
      prod(dpois(sim_Dx, lambda = sim_Ex*q_x))
    }


  }else{

    y = log(data$qx/(1-data$qx))
    sigma = sqrt(median(fit$post.samples$sigma2))

    likelihood <- function(par){
      gomp = par[1]*exp(par[2]*data$x)
      prod(dnorm(y, mean = log(gomp), sd = sigma))
    }

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
