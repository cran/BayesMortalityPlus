## Funçoes para calcular a curva, dado uma idade e os parâmetros dos modelos
hp_curve_9 = function(x, p){
  a = p[1] ; b = p[2] ; c = p[3] ; d = p[4] ; e = p[5] ; f = p[6] ; g = p[7] ; h = p[8] ; k = p[9]
  a^((x+b)^c) + d*exp(-e*log(x/f)^2) + g*h^x / (1 + k*g*h^x)
}

hp_curve <- function(x, p){
  a = p[1]; b = p[2]; c = p[3]; d = p[4]; e = p[5]; f = p[6]; g = p[7]; h = p[8]
  a ^ ((x + b) ^ c) + d * exp(-e * (log(x) - log(f))^2) + g * h ^ x
}

## Pequena função para verificar o número de decimais necessários na escala y
decimal <- function(n){
  cont <- 0
  while(n%%10 < 1){
    n <- 10*n; cont <- cont + 1
  }
  return(cont)
}

################# HP OPTIMIZATION
###### Funcoes de perda
### Curva HP de 9 parametros
f.perda9 = function(p, Ex, Dx, x){ ## e a funcao perda default -> O K e sempre positivo!!!
  mu = hp_curve_9(x, exp(p)); mx = Dx/Ex
  loss = log(mu/mx)^2
  loss[is.infinite(loss)] <- 10^5
  return(sum(loss, na.rm = T))
}

#### Essa versao permite estimar K nos negativos tambem
f.perda9.alt = function(p, Ex, Dx, x){
  p <- c(exp(p[-9]), p[9])
  mu = hp_curve_9(x, p); mx = Dx/Ex
  loss = log(mu/mx)^2
  loss[is.infinite(loss)] <- 10^5
  return(sum(loss, na.rm = T))
}

### Curva HP de 8 parametros
f.perda = function(p, Ex, Dx, x){
  mu = hp_curve(x, exp(p)); mx = Dx/Ex
  mu = mu/(1+mu)
  loss = log(mu/mx)^2
  loss[is.infinite(loss)] <- 10^5
  return(sum(loss, na.rm = T))
}


#### Funcao de otimizacao
optim_HP <- function(x, Ex, Dx, curve = c("8par", "9par")){
  start = c(5e-04, 0.004, 0.08, 0.001, 10, 17, 5e-05, 1.1, 1) ## Start do mortalitlyLaws
  if(curve == "8par"){
    start = start[-9]
    aux <- nlminb(start = log(start),
                  objective = f.perda, x = x, Ex = Ex, Dx = Dx,
                  lower = log(c(1e-16, 1e-16, 1e-16, 1e-16, 1e-16, 15+1e-16, 1e-16, 1)),
                  upper = log(c(1-1e-16, 1-1e-16, 1-1e-16, 1-1e-16, 100, 110-1e-16, 1-1e-16, 100)),
                  control = list(eval.max = 5000, iter.max = 5000))
    return(exp(aux$par))
  }else if(curve == "9par"){
    alt = FALSE
    if(alt){
      aux <- nlminb(start = log(start),
                    objective = f.perda9.alt, x = x, Ex = Ex, Dx = Dx,
                    lower = c(log(c(1e-16, 1e-16, 1e-16, 1e-16, 1e-16, 15+1e-16, 1e-16, 1)), -4),
                    upper = c(log(c(1-1e-16, 1-1e-16, 1-1e-16, 1-1e-16, 100, 110-1e-16, 1-1e-16, 100)), 100),
                    control = list(eval.max = 5000, iter.max = 5000))
      return(c(aux$par[-9], aux$par[9]))
    }else{
      aux <- nlminb(start = log(start),
                    objective = f.perda9, x = x, Ex = Ex, Dx = Dx,
                    lower = log(c(1e-16, 1e-16, 1e-16, 1e-16, 1e-16, 15+1e-16, 1e-16, 1, 1e-16)),
                    upper = log(c(1-1e-16, 1-1e-16, 1-1e-16, 1-1e-16, 100, 110-1e-16, 1-1e-16, 100, 100)),
                    control = list(eval.max = 5000, iter.max = 5000))
      return(exp(aux$par))
    }
  }
}
