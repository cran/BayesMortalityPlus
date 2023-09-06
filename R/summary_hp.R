#' @name summary.HP
#' @rdname summary.HP
#'
#' @title HP: Summary
#'
#' @description Summarizes information from the parameters' markov chains of a fitted `HP` or `ClosedHP` model.
#'
#'
#' @param object A `HP` or `ClosedHP` object, result of a call to hp() or hp_close() function.
#' @param ... Further arguments passed to or from other methods.
#'
#' @return A data.frame object with the mean, standard deviation and 2.5%, 50% and 97.5% quantiles of a fitted `HP` or `ClosedHP` model.
#'
#' @examples
#'  ## Importing mortality data from the USA available on the Human Mortality Database (HMD):
#' data(USA)
#'
#' ## Selecting the exposure and death count of the 2010 male population ranging from 0 to 90 years old
#' USA2010 = USA[USA$Year == 2010,]
#' x = 0:90
#' Ex = USA2010$Ex.Male[x+1]
#' Dx = USA2010$Dx.Male[x+1]
#'
#' ## Fitting binomial model
#' fit = hp(x = x, Ex = Ex, Dx = Dx, M = 5000, bn = 0, thin = 10)
#' summary(fit)
#'
#' @seealso [summary.DLM()] for `DLM` method.
#'
#' @export
summary.HP <- function(object, ...){ fit = object; return(fit$summary) }

#' @export
summary.ClosedHP <- function(object, ...){
  fit = object
  resumo = data.frame(
    x = min(fit$data$x):max(fit$data$x),
    x1 = apply(fit$qx, 2, mean),
    x2 = apply(fit$qx, 2, sd),
    t(apply(fit$qx, 2, quantile, probs = c(0.025, 0.5, 0.975))))
  colnames(resumo) <- c("age", "mean", "sd", "2.5%", "50.0%", "97.5%")
  return(resumo)
}
