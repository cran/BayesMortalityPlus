#' @name summary.DLM
#' @rdname summary.DLM
#'
#' @title Summary for DLM fitted models
#'
#' @description Summarizes information from the parameters' markov chains of a fitted `DLM` or `ClosedDLM` model.
#'
#'
#' @param object A `DLM` or `ClosedDLM` object, result of a call to dlm() or dlm_close() function.
#' @param digits An integer indicating the number of decimals places.
#' @param ... Further arguments passed to or from other methods.
#'
#' @return A data.frame object with the mean, standard deviation and 2.5%, 50% and 97.5% quantiles of a fitted `DLM` or `ClosedDLM` model.
#'
#' @examples
#' ## Importing mortality data from the USA available on the Human Mortality Database (HMD):
#' data(USA)
#'
#' ## Selecting the log mortality rate of the 2010 male population ranging from 0 to 100 years old
#' USA2010 = USA[USA$Year == 2010,]
#' x = 0:100
#' Ex = USA2010$Ex.Male[x+1]
#' Dx = USA2010$Dx.Male[x+1]
#' qx_t = Dx/Ex
#' qx_t = 1 - exp(-qx_t)
#' y = log(qx_t)
#'
#' ## Fitting DLM
#' fit = dlm(y, M = 100, bn = 20, thin = 1)
#' summary(fit)
#'
#' @seealso [summary.HP()] for `HP` method.
#'
#' @export
summary.DLM <- function(object, digits = 5, ...){
  fit = object
  summ = rbind(resumo(fit$sig2),
               resumo(fit$mu))
  # summ = as.data.frame(summary)
  row.names(summ) = c("sigma2", paste0("mu[", fit$info$ages, "]"))
  return(round(summ, digits))
}

#' @export
summary.ClosedDLM <- function(object, ...){
  fit = object
  summ = as.data.frame(cbind(min(fit$info$ages):max(fit$info$ages), resumo(fit$qx)))
  colnames(summ) <- c("age", "mean", "sd", "2.5%", "50.0%", "97.5%")
  return(summ)
}

#'
resumo <- function(x){
  if(is.matrix(x)){
    resumo = data.frame(
      x1 = apply(x, 2, mean),
      x2 = apply(x, 2, sd),
      t(apply(x, 2, quantile, probs = c(0.025, 0.5, 0.975))))
  }else{
    resumo = data.frame(
      x1 = mean(x),
      x2 = sd(x),
      t(quantile(x, probs = c(0.025, 0.5, 0.975))))
  }
  colnames(resumo) <- c("mean", "sd", "2.5%", "50.0%", "97.5%")
  return(resumo)
}
