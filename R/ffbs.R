## Filtering
ff = function(m0, C0, y,  alpha0, beta0, Ft, Gt, delta){

  N = length(y)
  p = length(m0)
  resultado.m = matrix(NA_real_, N, p)
  resultado.C = array(NA_real_, c(N,p,p))
  resultado.W = array(NA_real_, c(N,p,p))
  resultado.a = matrix(NA_real_, N, p)
  resultado.R = array(NA_real_, c(N,p,p))
  resultado.f = c()
  resultado.Q = c()
  resultado.alpha = c()
  resultado.beta = c()
  if(length(delta) == 1){ delta = rep(delta, N) }

  V = 0.01

  ## Kalman Filter
  ### Step 1

  Wt = C0 * (1 - delta[1]) / delta[1]

  at = Gt %*% m0
  Rt = Gt %*% C0 %*% t(Gt) + Wt
  ft = Ft %*% at
  Qt = (Ft %*% Rt %*% t(Ft) + V)[1,1]
  et = y[1] - ft
  At = (Rt %*% t(Ft)) / Qt
  mt = at + At %*% et
  Ct = Rt - At %*% Ft %*% Rt
  alphat = alpha0 + 1/2
  betat = beta0 + 0.5*t(et)%*%et/Qt


  resultado.m[1,] = mt
  resultado.C[1,,] = Ct
  resultado.W[1,,] = Ct * (1 - delta[1]) / delta[1]
  resultado.a[1,] = at
  resultado.R[1,,] = Rt
  resultado.f[1] = ft ##
  resultado.Q[1] = Qt ##
  resultado.alpha[1] = alphat
  resultado.beta[1] = betat


  ### Step 2
  for (j in 2:N) {

    Wt = Ct * (1 - delta[j]) / delta[j]
    at = Gt %*% mt
    Rt = Gt %*% Ct %*% t(Gt) + Wt
    ft = Ft %*% at
    Qt = (Ft %*% Rt %*% t(Ft) + V)[1,1]
    et = y[j] - ft
    At = (Rt %*% t(Ft)) / Qt
    mt = at + At %*% et
    Ct = Rt - At %*% Ft %*% Rt
    alphat = alphat + 1/2
    betat = betat + 0.5*t(et)%*%et/Qt

    resultado.m[j,] = mt
    resultado.C[j,,] = Ct
    resultado.W[j,,] = Wt
    resultado.a[j,] = at
    resultado.R[j,,] = Rt
    resultado.f[j] = ft ##
    resultado.Q[j] = Qt ##
    resultado.alpha[j] = alphat
    resultado.beta[j] = betat

  }

  return(list(m = resultado.m, C = resultado.C,
              a = resultado.a, R = resultado.R,
              W = resultado.W, f = resultado.f,
              Qt = resultado.Q,
              alpha = resultado.alpha,
              beta = resultado.beta))
}

## Backward Sampling
bs = function(m, C, a, R, Gt, alpha, beta){

  N = nrow(m)
  p = ncol(m)

  as = matrix(NA, N, p)
  Rs = array(NA, c(N, p, p))
  # theta <- matrix(NA,N,p)

  ## Distribuicao de thetaT
  as[N,] = m[N,] ##
  Rs[N,,] = C[N,,] ##
  alpha = alpha[N]
  beta = beta[N]

  ### step 3 - Smoothing
  for (t in (N - 1):1) {

    Bt = C[t,,] %*% t(Gt) %*% chol2inv(chol(R[t + 1,,]))
    Rs[t,,] = C[t,,] + Bt %*% (Rs[t + 1,,] - R[t + 1,,]) %*% t(Bt)
    as[t,] = m[t,] + Bt %*% (as[t + 1,] - a[t + 1,])
  }

  return(list(as = as, Rs = Rs, alpha = alpha, beta = beta))
}

# Filtering and smoothing with FFBS and discount factor W
ffbs <- function(m0, C0, y,  alpha0, beta0, Ft, Gt, delta){

  aux.f = ff(m0, C0, y,  alpha0, beta0, Ft, Gt, delta)

  res = bs(aux.f$m, aux.f$C, aux.f$a, aux.f$R, Gt, aux.f$alpha, aux.f$beta)
  # res$W = aux.f$W

  return(res)

}
