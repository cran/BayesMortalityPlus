# Filtragem com V fixo e fator de desconto para W
#
# @param m0 vetor de dimensao p com os valores iniciais de mt da filtragem
# @param C0 matriz de dimensao pxp com os valores iniciais de Ct da filtragem
# @param y vetor de tamanho N
# @param V variância observacional - escalar
# @param Ft0 matriz de dimensao Nxp
# @param Gt matriz de dimensao pxp
# @param delta fator de desconto - escalar
#
# @return Uma lista contendo mt, Ct, at, Rt e Wt
#
# @examples
# n = length(Nile)
# mld = ff(m0 = 0, C0 = 100000, y = Nile, V = 15099.8,
#          Ft0 = matrix(rep(1,n),n,1), Gt = 1, delta = 0.85)
# xg = 1:length(Nile)
# q = qnorm(1 - 0.05/2)
# qfinf = mld$m - q*sqrt(as.vector(mld$C))
# qfsup = mld$m + q*sqrt(as.vector(mld$C))
# ts.plot(as.vector(Nile))
# polygon(x = c(xg, rev(xg)), y = c(qfinf, rev(qfsup)),
#         col = adjustcolor('lightcoral', alpha.f = 0.3), border = NA)
# lines(as.vector(mld$m), col = 2, lwd = 2)
ff = function(m0, C0, y, V, Ft, Gt, delta){

  # if (!is.matrix(C0)) stop("C0 deve ser uma matriz")
  # if (dim(C0)[1] != dim(C0)[2]) stop("C0 deve ser uma matriz quadrada")
  # if (!is.vector(m0)) stop("m0 deve ser um vetor")
  # if (length(m0) != dim(C0)[1]) stop("Dimensões de m0 e C0 incompativeis")
  # if (!is.vector(y)) stop("y deve ser um vetor")
  # if (!is.matrix(Ft0)) stop("Ft0 deve ser uma matriz")
  # if (dim(Ft0)[1] != length(y)) stop("Dimensões de Ft0 e y imcompativeis")
  # if (!is.matrix(Gt)) stop("Gt deve ser uma matriz")
  # if (delta < 0 | delta > 1) stop("delta deve estar entre 0 e 1")

  N = nrow(y)
  p = length(m0)
  resultado.m = matrix(NA, N, p)
  resultado.C = array(NA,c(N,p,p))
  resultado.W = array(NA,c(N,p,p))
  resultado.a = matrix(NA, N, p)
  resultado.R = array(NA,c(N,p,p))

  ## Filtro de Kalman
  ### passo 1 (inicializacao)

  Wt = C0 * (1 - delta) / delta

  #  if(is.matrix(Ft0) == TRUE){ Ft = Ft0[1,] }
  at = Gt %*% m0
  Rt = Gt %*% C0 %*% t(Gt) + Wt
  ft = Ft %*% at
  Qt = Ft %*% Rt %*% t(Ft) + V
  et = y[1,] - ft
  At = Rt %*% t(Ft) %*% solve(Qt)
  mt = at + At %*% et   ### first moment
  Ct = Rt - At %*% Ft %*% Rt   ## second moment

  resultado.m[1,] = mt
  resultado.C[1,,] = Ct
  resultado.W[1,,] = Ct * (1 - delta) / delta
  resultado.a[1,] = at
  resultado.R[1,,] = Rt

  ### passo 2 (atualizacao)
  for (j in 2:N) {
    Wt = Ct * (1 - delta) / delta
    # if(is.matrix(Ft0) == TRUE){ Ft = Ft0[j,] }
    at = Gt %*% mt
    Rt = Gt %*% Ct %*% t(Gt) + Wt
    ft = Ft %*% at
    Qt = Ft %*% Rt %*% t(Ft) + V
    et = y[j,] - ft
    At = Rt %*% t(Ft) %*% solve(Qt)
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

# Backward Sampling
#
# Usada internamente na funcao ffbs()
#
# @param m vetor de dimensao p com o resultado da filtragem
# @param C matriz de dimensao pxp com o resultado da filtragem
# @param a vetor de dimensao p com o resultado da filtragem
# @param R matriz de dimensao pxp com o resultado da filtragem
# @param Gt matriz de dimensao pxp
#
# @return Uma lista contendo theta, mt, Ct, as e Rs
#
# @importFrom MASS mvrnorm
#
bs = function(m,C,a,R,Gt){

  N = nrow(m)
  p = ncol(m)

  as = matrix(NA, N, p)
  Rs = array(NA,c(N,p,p))
  theta <- matrix(NA,N,p)

  as[N,] = m[N,]
  Rs[N,,] = C[N,,]

  ### draw theta_T - page 162 petris petroni
  theta[N,] <- MASS::mvrnorm(1, as[N,], Rs[N,,])

  ### step 3 - algorithm 4.1 Backward Sampling
  for (t in (N - 1):1) {

    Bt = C[t,,] %*% t(Gt) %*% solve(R[t + 1,,])

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

# Filtragem e suavização com FFBS e fator de desconto para W
#
# Usada internamente na funcao gibbsSigma2
#
# @param m0 vetor de dimensao p com os valores iniciais de mt da filtragem
# @param C0 matriz de dimensao pxp com os valores iniciais de Ct da filtragem
# @param y vetor de tamanho N
# @param V variância observacional - escalar
# @param Ft0 matriz de dimensao Nxp
# @param Gt matriz de dimensao pxp
# @param delta fator de desconto - escalar
#
# @return Uma lista contendo theta, mt, Ct, as, Rs e W
#
ffbs <- function(m0, C0, y, V, Ft, Gt, delta){

  aux.f = ff(m0, C0, y, V, Ft, Gt, delta)

  res = bs(aux.f$m,aux.f$C,aux.f$a,aux.f$R,Gt)
  res$W = aux.f$W

  return(res)

}


# Estimação de V constante via Gibss
#
# Utlizada ffbs e um passo de Gibbs para estimar uma variancia constante no tmepo.
#
# @param m0 vetor de dimensao p com os valores iniciais de mt da filtragem
# @param C0 matriz de dimensao pxp com os valores iniciais de Ct da filtragem
# @param y vetor de tamanho N
# @param Ft0 matriz de dimensao Nxp
# @param Gt matriz de dimensao pxp
# @param delta fator de desconto - escalar
# @param sig2k valor inicial para sigma2
# @param alpha priori inversa gama para variancia (ou precisao, confirmar)
# @param beta priori inversa gama para variancia (ou precisao, confirmar)
# @param nit numero de iteracoes
# @param shiny habilita as funcoes para exibir progresso no shiny e cancelar o calculo
# @param status_file só é usado quando shiny = T.
#
# @return Uma lista contendo uma cadeia da posteriori de mu, beta e sigma2.
#
# @examples
#
# data_hp$qx <- data_hp$dx/data_hp$nx
# data_hp$qx <- 1 - exp(-data_hp$qx)
# y <- log(data_hp$qx)
# plot(y, t='l')
# n <- length(y)
#
# Gt <- matrix(c(1,1,0,1),  ncol=2, byrow=TRUE)
# Ft0 <- matrix(c(1,0), n , ncol=2, byrow=TRUE)
# m0 <-  rep(0, 2)
# C0 <- diag(100,2)
# delta <- 0.65
# V <- 0.002674713
#
# res <- gibbsSigma2(m0,C0,y,Ft0,Gt,delta,0.002674713,0.01,0.01,500)
#
# plot(res$sig2, t='l')
#
# mu <- apply(res$mu, 2, median)
# Ft <- t(c(1,0))
# media <- mu
# med = exp(media)[2:76]
# plot(data_hp$x, med, type = "l", log = "y", ylim=c(8e-05, 1e-01))
# points(data_hp$x,data_hp$qx,pch = 19, cex = 0.5, col = "red")
#
# @import progress
#
gibbsSigma2 <- function(m0,C0,y,Ft0,Gt,delta,sig2k,alpha,beta,nit,shiny = F, status_file = NULL){
  n <- length(y)
  p = length(m0)
  Ft <- t(as.matrix(Ft0[1,]))
  sig2.post <- NULL
  mu.post <- matrix(NA, nit, ncol = n)
  # theta.post <- matrix(NA, nit, ncol = n)
  theta.post <- array(dim = c(nit, n, ncol(Ft0)))
  Wt = array(NA,dim = c(n,p,p,nit))

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
