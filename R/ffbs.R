## Filtering
ff = function(m0, C0, y, V, Ft, Gt, delta){

  N = nrow(y)
  p = length(m0)
  resultado.m = matrix(NA_real_, N, p)
  resultado.C = array(NA_real_,c(N,p,p))
  resultado.W = array(NA_real_,c(N,p,p))
  resultado.a = matrix(NA_real_, N, p)
  resultado.R = array(NA_real_,c(N,p,p))

  ## Kalman Filter
  ### Step 1

  Wt = C0 * (1 - delta) / delta

  #  if(is.matrix(Ft0) == TRUE){ Ft = Ft0[1,] }
  at = Gt %*% m0
  Rt = Gt %*% C0 %*% t(Gt) + Wt
  ft = Ft %*% at
  Qt = Ft %*% Rt %*% t(Ft) + V
  et = y[1,] - ft
  At = Rt %*% t(Ft) %*% chol2inv(chol(Qt))
  mt = at + At %*% et   ### first moment
  Ct = Rt - At %*% Ft %*% Rt   ## second moment

  resultado.m[1,] = mt
  resultado.C[1,,] = Ct
  resultado.W[1,,] = Ct * (1 - delta) / delta
  resultado.a[1,] = at
  resultado.R[1,,] = Rt

  ### Step 2
  for (j in 2:N) {
    Wt = Ct * (1 - delta) / delta
    # if(is.matrix(Ft0) == TRUE){ Ft = Ft0[j,] }
    at = Gt %*% mt
    Rt = Gt %*% Ct %*% t(Gt) + Wt
    ft = Ft %*% at
    Qt = Ft %*% Rt %*% t(Ft) + V
    et = y[j,] - ft
    At = Rt %*% t(Ft) %*% chol2inv(chol(Qt))
    mt = at + At %*% et  ### mean
    Ct = Rt - At %*% Ft %*% Rt ### variance

    resultado.m[j,] = mt
    resultado.C[j,,] = Ct
    resultado.W[j,,] = Wt
    resultado.a[j,] = at
    resultado.R[j,,] = Rt
  }

  return(list(m = resultado.m, C = resultado.C, a = resultado.a, R = resultado.R, W = resultado.W))
}

## Backward Sampling
bs = function(m,C,a,R,Gt){

  N = nrow(m)
  p = ncol(m)

  as = matrix(NA_real_, N, p)
  Rs = array(NA_real_,c(N,p,p))
  theta <- matrix(NA_real_,N,p)

  as[N,] = m[N,]
  Rs[N,,] = C[N,,]

  ### draw theta_T - page 162 petris petroni
  theta[N,] <- MASS::mvrnorm(1, as[N,], Rs[N,,])

  ### step 3 - algorithm 4.1 Backward Sampling
  for (t in (N - 1):1) {

    Bt = C[t,,] %*% t(Gt) %*% chol2inv(chol(R[t + 1,,]))

    # Rs[t,,] = C[t,,] + Bt %*% (Rs[t + 1,,] - R[t + 1,,]) %*% t(Bt)

    # as[t,] = m[t,] + Bt %*% (as[t + 1,] - a[t + 1,])

    ht <- m[t,] + Bt %*% (theta[t + 1,] - a[t + 1,])
    Ht = C[t,,] - Bt %*% R[t + 1,,] %*% t(Bt)
    ### draw theta_t
    theta[t,]  = MASS::mvrnorm(1,ht, Ht)
    # theta[t,]  = MASS::mvrnorm(1,as[t,], Rs[t,,])
  }
  return(list(d = theta, m = m, C = C))
  #m.m = as, C.C = Rs))

}

# Filtering and smoothing with FFBS and discount factor W
ffbs <- function(m0, C0, y, V, Ft, Gt, delta){

  aux.f = ff(m0, C0, y, V, Ft, Gt, delta)

  res = bs(aux.f$m,aux.f$C,aux.f$a,aux.f$R,Gt)
  res$W = aux.f$W

  return(res)

}


# Sampling V via Gibbs
gibbsSigma2 <- function(m0,C0,y,Ft0,Gt,delta,sig2k,alpha,beta,nit,shiny = F, status_file = NULL){
  n <- length(y)
  p = length(m0)
  Ft <- t(as.matrix(Ft0[1,]))
  sig2.post <- NULL
  mu.post <- matrix(NA_real_, nit, ncol = n)
  # theta.post <- matrix(NA, nit, ncol = n)
  theta.post <- array(dim = c(nit, n, ncol(Ft0)))
  Wt = array(NA_real_,dim = c(n,p,p,nit))

  pb  = progress::progress_bar$new(format = "Simulating [:bar] :percent in :elapsed",
                                   total = nit, clear = FALSE, width = 60)

  for (k in 1:nit) {

    pb$tick()

    V = sig2k

    ## FFBS for thetas (each age x)
    mld <- ffbs(m0 , C0 , y, V = V, Ft0 , Gt , delta)
    dt <- mld$d
    Wt[,,,k] = mld$W

    mu.post[k,] <- Ft%*%t(dt) ;   theta.post[k,,] <- dt ;
    muk <- mu.post[k,]

    ####### gibbs for sigma2
    alpha.star <- alpha + (n/2)
    beta.star <- 0.5*sum((y-muk)^2) + beta ;
    sig2.post[k] <- 1/rgamma(1, alpha.star, beta.star)
    sig2k <- sig2.post[k]
  }
  return(list(mu = mu.post, theta = theta.post, sig2 = sig2.post, Wt = Wt))
}
