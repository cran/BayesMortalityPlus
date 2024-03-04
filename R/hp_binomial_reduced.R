hp_binomial_red <- function(x, Ex, Dx, M = 100000, bn = round(M/4), thin = 1, m = NULL, v = NULL,
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
  a = b = c = 0 #reduced model
  inits[1:3] = 0

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
    logvero = sum(Dx*log(qx), na.rm = T) + sum((Ex - Dx)*log(1 - qx), na.rm = T);
    return(logvero)
  }

  ### Jacobian
  jac_logit = function(x) - log(abs(x - x^2))
  jac_log = function(x) - log(x)
  jac_f = function(x){
    - log(110 - x) - log(x - 15)
  }

  ##############################################################################################

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

  ### proposed parameters for the block (d, e, f, g, h)
  U = diag(5)
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

    ##### 'd', 'e', 'f', 'g', 'h' estimated in one block (joint)

    ### Covariance matrix (Metropolis-Hastings)
    if(k < 1000){
      V = sd*U
    }else if(k%%10 == 0) { ### updating every 10 iterations
      X = theta.post[c((k-1000):(k-1)), 4:8]
      X[,1] = log(X[,1]/(1 - X[,1]))                                      ## D
      X[,2] = log(X[,2])                                                  ## E
      X[,3] = (X[,3] - 15)/(110 - 15); X[,3] = log(X[,3]/(1 - X[,3]))     ## F
      X[,4] = log(X[,4]/(1 - X[,4]))                                      ## G
      X[,5] = log(X[,5])                                                  ## H
      V = sd*(eps*U + var(X))
    }

    aux = theta.post[k-1, 4:8]
    aux[1] = log(aux[1]/(1 - aux[1]))                                      ## D
    aux[2] = log(aux[2])                                                   ## E
    aux[3] = (aux[3] - 15)/(110 - 15); aux[3] = log(aux[3]/(1 - aux[3]))   ## F
    aux[4] = log(aux[4]/(1 - aux[4]))                                      ## G
    aux[5] = log(aux[5])                                                   ## H
    prop = MASS::mvrnorm(1, mu = aux, Sigma = V)
    d.prop = exp(prop[1]) / (1 + exp(prop[1]))
    e.prop = exp(prop[2])
    f.prop = 15 + (110 - 15)*(exp(prop[3])/ (1 + exp(prop[3])))
    g.prop = exp(prop[4]) / (1 + exp(prop[4]))
    h.prop = exp(prop[5])

    lverok    = like.HPBayes(Ex, Dx, x, theta.post[k-1,1], theta.post[k-1,2], theta.post[k-1,3], theta.post[k-1,4], theta.post[k-1,5], theta.post[k-1,6], theta.post[k-1,7], theta.post[k-1,8], theta.post[k-1,9])
    lveroprop = like.HPBayes(Ex, Dx, x, a                , b                , c                , d.prop           , e.prop           , f.prop           , g.prop           , h.prop           , theta.post[k-1,9])


    auxk    = lverok    + dbeta(theta.post[k-1,4], alpha.d, beta.d, log = T) + dgamma(theta.post[k-1,5], alpha.e, beta.e, log = T) + dnorm(theta.post[k-1,6], media.f, sqrt(variancia.f), log = T) + dbeta(theta.post[k-1,7], alpha.g, beta.g, log = T) + dgamma(theta.post[k-1,8], alpha.h, beta.h, log = T) + jac_logit(theta.post[k-1,4]) + jac_log(theta.post[k-1,5]) + jac_f(theta.post[k-1,6]) + jac_logit(theta.post[k-1,7]) + jac_log(theta.post[k-1,8])
    auxprop = lveroprop + dbeta(d.prop           , alpha.d, beta.d, log = T) + dgamma(e.prop           , alpha.e, beta.e, log = T) + dnorm(f.prop           , media.f, sqrt(variancia.f), log = T) + dbeta(g.prop           , alpha.g, beta.g, log = T) + dgamma(h.prop           , alpha.h, beta.h, log = T) + jac_logit(d.prop)            + jac_log(e.prop)            + jac_f(f.prop)            + jac_logit(g.prop)            + jac_log(h.prop)

    ratio = auxprop - auxk; test = runif(1)
    if(is.na(ratio) || is.nan(ratio)){ ratio = -Inf }
    if (ratio > log(test)) {
      theta.post[k,4] = d.prop
      theta.post[k,5] = e.prop
      theta.post[k,6] = f.prop
      theta.post[k,7] = g.prop
      theta.post[k,8] = h.prop
      cont[4:8] = cont[4:8] + 1
    } else {
      theta.post[k,4:8] = theta.post[k-1,4:8]
    }

    theta.post[k, 1] = a
    theta.post[k, 2] = b
    theta.post[k, 3] = c
    theta.post[k,9] = theta.post[k-1,9] ### fixed K

    limite_inf = which(theta.post[k,4:8] <= c(1e-16, 1e-16, 15+1e-16, 1e-16, 1e-16))
    limite_sup = which(theta.post[k,4:8] >= c(1-1e-5, Inf, 110-1e-5, 1-1e-5, Inf))

    if(length(limite_inf) > 0){
      theta.post[k, 3 + limite_inf] <- theta.post[k-1, 3 + limite_inf]
      if(k > bn){
        param_problemas = append(param_problemas, param[3 + limite_inf]); warn = T
      }
    }

    if(length(limite_sup) > 0){
      theta.post[k, 3 + limite_sup] <- theta.post[k-1, 3 + limite_sup]
      if(k > bn){
        param_problemas = append(param_problemas, param[3 + limite_sup]); warn = T
      }
    }


  })

  if(warn){ warning(paste0("MCMC may have had some issues with the parameter(s): ", paste(sort(unique(param_problemas)), collapse = ", "), ".\nCheck the 'plot_chain' function output to visualize the parameters chain. It might be helpful to assign informative prior distribution for these parameters. See ?hp.")) }

  ##############################################################################################
  ### Return

  ## Final samples
  theta.post = theta.post[seq(bn+1+1, M+1, by = thin),]

  return(list(theta.post = theta.post, sigma2 = NULL, cont = cont, inits = inits))

}
