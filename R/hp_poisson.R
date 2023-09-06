hp_poisson <- function(x, Ex, Dx, M = 100000, bn = round(M/4), thin = 1, m = NULL, v = NULL,
                       inits = NULL, prop.control = NULL, K = NULL){

  ##############################################################################################
  #### optimal value for the inits
  Mx = Dx/Ex
  opt = optim_HP(x = x, Ex, Dx, curve = "9par")

  ## init validation
  if(is.null(inits)){

    inits[1] = opt[1]; inits[2] = opt[2]; inits[3] = opt[3]; inits[4] = opt[4]
    inits[5] = opt[5]; inits[6] = opt[6]; inits[7] = opt[7]; inits[8] = opt[8]

    if(inits[1] <= 0 || inits[1] >= 1){ inits[1] = 0.01 }
    if(inits[2] <= 0 || inits[2] >= 1){ inits[2] = 0.1 }
    if(inits[3] <= 0 || inits[3] >= 1){ inits[3] = 0.1 }
    if(inits[4] <= 0 || inits[4] >= 1){ inits[4] = 0.01 }
    if(inits[5] <= 0){ inits[5] = 5 }
    if(inits[6] <= 15 || inits[6] >= 110){ inits[6] = 25 }
    if(inits[7] <= 0 || inits[7] >= 1){ inits[7] = 0.01 }
    if(inits[8] <= 0){ inits[8] = 1.1 }

  }else if(length(inits) != 8 || any(is.na(inits[1:8]))){
    stop("inits vector with missing values (NA) or length not equal to 8.")
  }else if(any(inits[c(1:5,7,8)] <= 0) || any(inits[c(1:4,7)] >= 1) || inits[6] <= 15 || inits[6] >= 110){
    stop("Initial value(s) outside the range of possible values for the paramater(s).")
  }

  ## K
  k2 = ifelse(is.null(K), opt[9], K)

  ##############################################################################################
  ### auxs

  mu = function(x, a, b, c, d, e, f, g, h, k){
    ### HP function for specific age and parameters
    media = a^((x+b)^c) + d*exp(-e*(log(x)-log(f))^2) + (g*h^x)/(1 + k*g*h^x)
    return(media)
  }

  ### Log-likelihood (MCMC)
  like.HPBayes= function(Ex, Dx, x, a,b,c,d,e,f,g,h,k){

    qx = mu(x, a, b, c, d, e, f, g, h, k)
    # qx = 1 - exp(-mx)
    qx[qx < 1e-16] = 1e-16            # avoid numeric error w log(qx)
    qx[qx > 1 - 1e-16] = 1 - 1e-16    # avoid numeric error w log(1 - qx)
    # Ex = trunc(Ex)
    # logvero = sum(Dx*log(qx), na.rm = T) + sum((Ex - Dx)*log(1 - qx), na.rm = T)
    logvero = sum(Dx*log(Ex*qx) - Ex*qx, na.rm = T);
    return(logvero)
  }

  ### Jacobian
  jac_logit = function(x) - log(abs(x - x^2))
  jac_log = function(x) - log(x)
  jac_f = function(x){
    - log(110 - x) - log(x - 15)
  }

  ##############################################################################################

  ## Prioris
  if(v[1] < m[1]*(1-m[1])){
    alpha.a = ((1 - m[1])/v[1] - 1/m[1])*m[1]^2
    beta.a = alpha.a*(1/m[1] - 1)
  }

  if(v[2] < m[2]*(1-m[2])){
    alpha.b = ((1 - m[2])/v[2] - 1/m[2])*m[2]^2
    beta.b = alpha.b*(1/m[2] - 1)
  }

  if(v[3] < m[3]*(1-m[3])){
    alpha.c = ((1 - m[3])/v[3] - 1/m[3])*m[3]^2
    beta.c = alpha.c*(1/m[3] - 1)
  }

  if(v[4] < m[4]*(1-m[4])){
    alpha.d = ((1 - m[4])/v[4] - 1/m[4])*m[4]^2
    beta.d = alpha.d*(1/m[4] - 1)
  }

  alpha.e = (m[5]^2)/v[5]
  beta.e = m[5]/v[5]

  media.f = m[6]; variancia.f = v[6]

  if(v[7] < m[7]*(1-m[7])){
    alpha.g = ((1 - m[7])/v[7] - 1/m[7])*m[7]^2
    beta.g = alpha.g*(1/m[7] - 1)
  }

  alpha.h = (m[8]^2)/v[8]
  beta.h = m[8]/v[8]

  ## SD for prop. distributions
  sd = 1*prop.control

  ### proposed parameters for the block (a, b, c, d, e, f, g, h)
  U = diag(8)
  eps = 1e-10

  ##############################################################################################
  ## aux objects
  param = c("A", "B", "C", "D", "E", "F", "G", "H", "K")
  param_problemas = NULL; warn = F

  theta.post = matrix(NA_real_, ncol = 9, nrow = M + 1)

  ### Acceptance rates
  cont = rep(0, 9)

  ### Initial values
  theta.post[1,] = c(inits, k2)

  ## progress bar
  pb  = progress::progress_bar$new(format = "Simulating [:bar] :percent in :elapsed",total = M, clear = FALSE, width = 60)

  ##############################################################################################
  ## Fits

  system.time(for (k in 2:(M+1)) {

    pb$tick()

    ##### ''a', 'b', 'c', 'd', 'e', 'f', 'g', 'h' estimated in one block (joint)

    ### Covariance matrix (Metropolis-Hastings)
    if(k < 1000){
      V = sd*U
    }else if(k%%10 == 0) { ### updating every 10 iterations
      X = theta.post[c((k-1000):(k-1)), 1:8]
      X[,1] = log(X[,1]/(1 - X[,1]))                                      ## A
      X[,2] = log(X[,2]/(1 - X[,2]))                                      ## B
      X[,3] = log(X[,3]/(1 - X[,3]))                                      ## C
      X[,4] = log(X[,4]/(1 - X[,4]))                                      ## D
      X[,5] = log(X[,5])                                                  ## E
      X[,6] = (X[,6] - 15)/(110 - 15); X[,6] = log(X[,6]/(1 - X[,6]))     ## F
      X[,7] = log(X[,7]/(1 - X[,7]))                                      ## G
      X[,8] = log(X[,8])                                                  ## H
      V = sd*(eps*U + var(X))
    }

    aux = theta.post[k-1, 1:8]
    aux[1] = log(aux[1]/(1 - aux[1]))                                      ## A
    aux[2] = log(aux[2]/(1 - aux[2]))                                      ## B
    aux[3] = log(aux[3]/(1 - aux[3]))                                      ## C
    aux[4] = log(aux[4]/(1 - aux[4]))                                      ## D
    aux[5] = log(aux[5])                                                   ## E
    aux[6] = (aux[6] - 15)/(110 - 15); aux[6] = log(aux[6]/(1 - aux[6]))   ## F
    aux[7] = log(aux[7]/(1 - aux[7]))                                      ## G
    aux[8] = log(aux[8])                                                   ## H
    prop = MASS::mvrnorm(1, mu = aux, Sigma = V)
    a.prop = exp(prop[1]) / (1 + exp(prop[1]))
    b.prop = exp(prop[2]) / (1 + exp(prop[2]))
    c.prop = exp(prop[3]) / (1 + exp(prop[3]))
    d.prop = exp(prop[4]) / (1 + exp(prop[4]))
    e.prop = exp(prop[5])
    f.prop = 15 + (110 - 15)*(exp(prop[6])/ (1 + exp(prop[6])))
    g.prop = exp(prop[7]) / (1 + exp(prop[7]))
    h.prop = exp(prop[8])

    lverok    = like.HPBayes(Ex, Dx, x, theta.post[k-1,1], theta.post[k-1,2], theta.post[k-1,3], theta.post[k-1,4], theta.post[k-1,5], theta.post[k-1,6], theta.post[k-1,7], theta.post[k-1,8], theta.post[k-1,9])
    lveroprop = like.HPBayes(Ex, Dx, x, a.prop           , b.prop           , c.prop           , d.prop           , e.prop           , f.prop           , g.prop           , h.prop           , theta.post[k-1,9])


    auxk    = lverok    + dbeta(theta.post[k-1,1], alpha.a, beta.a, log = T) + dbeta(theta.post[k-1,2], alpha.b, beta.b, log = T) + dbeta(theta.post[k-1,3], alpha.c, beta.c, log = T) + dbeta(theta.post[k-1,4], alpha.d, beta.d, log = T) + dgamma(theta.post[k-1,5], alpha.e, beta.e, log = T) + dnorm(theta.post[k-1,6], media.f, sqrt(variancia.f), log = T) + dbeta(theta.post[k-1,7], alpha.g, beta.g, log = T) + dgamma(theta.post[k-1,8], alpha.h, beta.h, log = T) + jac_logit(theta.post[k-1,1]) + jac_logit(theta.post[k-1,2]) + jac_logit(theta.post[k-1,3]) + jac_logit(theta.post[k-1,4]) + jac_log(theta.post[k-1,5]) + jac_f(theta.post[k-1,6]) + jac_logit(theta.post[k-1,7]) + jac_log(theta.post[k-1,8])
    auxprop = lveroprop + dbeta(a.prop           , alpha.a, beta.a, log = T) + dbeta(b.prop           , alpha.b, beta.b, log = T) + dbeta(c.prop           , alpha.c, beta.c, log = T) + dbeta(d.prop           , alpha.d, beta.d, log = T) + dgamma(e.prop           , alpha.e, beta.e, log = T) + dnorm(f.prop           , media.f, sqrt(variancia.f), log = T) + dbeta(g.prop           , alpha.g, beta.g, log = T) + dgamma(h.prop           , alpha.h, beta.h, log = T) + jac_logit(a.prop)            + jac_logit(b.prop)            + jac_logit(c.prop)            + jac_logit(d.prop)            + jac_log(e.prop)            + jac_f(f.prop)            + jac_logit(g.prop)            + jac_log(h.prop)

    ratio = auxprop - auxk; test = runif(1)
    if(is.na(ratio) || is.nan(ratio)){ ratio = -Inf }
    if (ratio > log(test)) {
      theta.post[k,1] = a.prop
      theta.post[k,2] = b.prop
      theta.post[k,3] = c.prop
      theta.post[k,4] = d.prop
      theta.post[k,5] = e.prop
      theta.post[k,6] = f.prop
      theta.post[k,7] = g.prop
      theta.post[k,8] = h.prop
      cont[1:8] = cont[1:8] + 1
    } else {
      theta.post[k,1:8] = theta.post[k-1,1:8]
    }

    theta.post[k,9] = theta.post[k-1,9]

    limite_inf = which(theta.post[k,] <= c(1e-16, 1e-16, 1e-16, 1e-16, 1e-16, 15+1e-16, 1e-16, 1e-16, -Inf))
    limite_sup = which(theta.post[k,] >= c(1-1e-5, 1-1e-5, 1-1e-5, 1-1e-5, Inf, 110-1e-5, 1-1e-5, Inf, Inf))

    if(length(limite_inf) > 0){
      theta.post[k, limite_inf] <- theta.post[k-1, limite_inf]
      param_problemas = append(param_problemas, param[limite_inf]); warn = T
    }

    if(length(limite_sup) > 0){
      theta.post[k, limite_sup] <- theta.post[k-1, limite_sup]
      param_problemas = append(param_problemas, param[limite_sup]); warn = T
    }

  })

  if(warn){ warning(paste0("MCMC had some problem with parameter(s): ", paste(sort(unique(param_problemas)), collapse = ", "), ". Try to assign informative prior distributions.")) }

  ##############################################################################################
  ### Return

  ## Final samples
  theta.post = theta.post[seq(bn+1+1, M+1, by = thin),]

  return(list(theta.post = theta.post, sigma2 = NULL, cont = cont, inits = inits))

}
