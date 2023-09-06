#' @name expectancy.HP
#' @rdname expectancy.HP
#'
#' @title HP: Life expectancy
#'
#' @description This function computes the life expectancy for each age for Heligman-Pollard model.
#'
#'
#' @param x Object of the class `HP` or `ClosedHP` fitted by hp() or hp_close() functions.
#' @param Ex Numeric vector with the exposure by age. This argument is only necessary when using poisson and binomial models with objects of the class `HP`.
#' @param age Numeric vector specifying the ages to calculate the life expectancy. The default is a sequence (0, 10, 20, ...) until the last decade used in the fitted model.
#' @param graph Logical value (TRUE ou FALSE). If TRUE, it returns a plot.
#' @param max_age Positive number indicating the last age to be considered to compute the life expectancy (extrapolation will be considered to match the age interval if needed). This argument is only necessary with objects of the class `HP`.
#' @param prob A percentage specifying the probability of credible interval.
#' @param ... Further arguments passed to or from other methods.
#'
#'
#' @return A data.frame and (if graph = TRUE) a plot.
#'
#' @examples
#' ## Importing mortality data from the USA available on the Human Mortality Database (HMD):
#' data(USA)
#'
#' # Example 1: --------------------------------
#'
#' USA1990 = USA[USA$Year == 1990,]
#'
#' Ex = USA1990$Ex.Total[1:91]
#' Dx = USA1990$Dx.Total[1:91]
#' x = 0:90
#'
#' fit <- hp(x, Ex, Dx, model = "binomial", M = 1000, bn = 0, thin = 10)
#' expectancy(fit)
#'
#'
#' # Example 2: -------------------------------
#'
#' # Using some arguments:
#'
#' Ex = USA1990$Ex.Total[1:106]
#'
#' expectancy(fit, Ex = Ex, age = c(0,20,30,60,105),
#' max_age = 105, prob = 0.99, graph = FALSE)
#'
#'
#' @include fitted_hp.R
#'
#' @import ggplot2
#'
#' @seealso [expectancy.DLM()] and [expectancy.BLC()] for `DLM` and `BLC`  methods.
#'
#' [Heatmap.HP()] and [Heatmap.list()] for `HP` or `list` methods to drawing a Heatmap for the truncated life expectancy.
#'
#' @export
expectancy.HP <- function(x, Ex = NULL, age = NULL, graph = TRUE,
                          max_age = 110, prob = 0.95, ...){
  fit = x
  if(is.null(age)){ age = seq(0, max(fit$data$x),by = 10) }
  ## Checking age
  if(max(age) > max_age){
    stop("Invalid age interval. Check the max_age argument")
  }

  ## calculating qx and px ----
  # extrapolating the age interval to max_age

  if(fit$info$model %in% c("binomial","poisson")){
    if(is.null(Ex)){
      Ex <- c(fit$data$Ex, rep(fit$data$Ex[length(fit$data$Ex)], (max_age+1)-length(fit$data$Ex)))

      aux <- fitted(fit, age = 0:max_age, Ex = Ex, prob = prob)
      est_IC <- aux
    }else{
      aux <- fitted(fit, age = 0:max_age, Ex = Ex, prob = prob)
      est_IC <- aux
    }
  }else{
    aux <- fitted(fit, age = 0:max_age, prob = prob)
    est_IC <- aux
  }
  qx_est <- aux$qx_fitted


  exp_total <- rep(NA_real_, max_age); exp_inf <- rep(NA_real_,max_age); exp_sup <- rep(NA_real_,max_age)

  # cumprod for life expectancy (px)
  for (i in 1:max_age){
    exp_total[i] <- sum(cumprod(1-qx_est[i:max_age])) ## px
    exp_sup[i] <- sum(cumprod(1-est_IC$qi[i:max_age])) ## upper CI
    exp_inf[i] <- sum(cumprod(1-est_IC$qs[i:max_age])) ## lower CI
  }
  exp_total <- round(exp_total,2)
  exp_sup <- round(exp_sup,2)
  exp_inf <- round(exp_inf,2)


  tab <- data.frame(x = 0:max(age),
                    exp_total[1:(max(age)+1)],
                    exp_inf[1:(max(age)+1)],
                    exp_sup[1:(max(age)+1)])
  tab[is.na(tab)] = 0
  colnames(tab) <- c("Age","Expectancy","Lower CI","Upper CI")

  if(graph == TRUE){
    p <-  ggplot(data=tab) + theme_light() +
      geom_line(aes(x=Age,y=Expectancy)) +
      geom_ribbon(aes(x=Age, ymin=`Lower CI`, ymax=`Upper CI`), alpha=0.3)
    return(list(expectancy=tab[tab$Age %in% age,],
                plot=p))
  }else{
    return(tab[tab$Age %in% age,])
  }
}


#' @export
#'
expectancy.ClosedHP <- function(x, age = seq(0, max(fit$data$x),by = 10),
                                graph = TRUE, prob = 0.95, ...){
  fit = x
  max_age <- max(fit$data$x)
  ###sanity
  if(max(age) > max_age){
    stop("Invalid age interval. Check the ages modeled")
  }

  ## calculating qx and px
  aux <- fitted(fit, prob = prob)
  qx_est <- aux$qx_fitted
  est_IC <- aux

  exp_total <- rep(NA_real_, max_age); exp_inf <- rep(NA_real_,max_age); exp_sup <- rep(NA_real_,max_age)

  # cumprod for life expectancy (px)
  for (i in 1:max_age){
    exp_total[i] <- sum(cumprod(1-qx_est[i:max_age])) ## px
    exp_sup[i] <- sum(cumprod(1-est_IC$qi[i:max_age])) ## upper CI
    exp_inf[i] <- sum(cumprod(1-est_IC$qs[i:max_age])) ## lower CI
  }
  exp_total <- round(exp_total,2)
  exp_sup <- round(exp_sup,2)
  exp_inf <- round(exp_inf,2)


  tab <- data.frame(x = 0:(max(age)),
                    exp_total[1:(max(age)+1)],
                    exp_inf[1:(max(age)+1)],
                    exp_sup[1:(max(age)+1)])
  tab[is.na(tab)] = 0
  colnames(tab) <- c("Age","Expectancy","Lower CI","Upper CI")

  if(graph == TRUE){
    p <-  ggplot(data=tab) + theme_light() +
      geom_line(aes(x=Age,y=Expectancy)) +
      geom_ribbon(aes(x=Age, ymin=`Lower CI`, ymax=`Upper CI`), alpha=0.3)
    return(list(expectancy=tab[tab$Age %in% age,],
                plot=p))
  }else{
    return(tab[tab$Age %in% age,])
  }
}
