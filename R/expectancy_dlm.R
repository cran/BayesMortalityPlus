#' @name expectancy.DLM
#' @rdname expectancy.DLM
#'
#' @title DLM: Life expectancy
#'
#' @description This function computes the life expectancy for each age for Dynamic Linear model.
#'
#'
#' @param x Object of the following classes: `DLM` or `ClosedDLM`.
#' @param age Numeric vector specifying the ages to calculate the life expectancy. The default is a sequence (0, 10, 20, ...) until the last decade used in the fitted model.
#' @param graph Logical value (TRUE ou FALSE). If TRUE, it also returns a plot. The default value is TRUE.
#' @param max_age Positive number indicating the last age to be considered to compute the life expectancy (prediction will be considered to match the age interval if needed). This argument is only necessary with objects of the class `DLM`.
#' @param prob A number specifying the probability of credible interval. The default value is 0.95.
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
#' Ex = USA1990$Ex.Total[1:111]
#' Dx = USA1990$Dx.Total[1:111]
#'
#' qx_t <- Dx/Ex
#' qx_t <- 1 - exp(-qx_t)
#' y <- log(qx_t)
#'
#' fit <- dlm(y, M = 100, bn = 20, thin = 1)
#' expectancy(fit)
#'
#' # Example 2: -------------------------------
#'
#' # Using some arguments:
#'
#' expectancy(fit, age = c(0,20,30,60),
#' prob = 0.99, max_age = 90, graph = FALSE)
#'
#'
#' @include qx_ci.R
#' @include fitted_dlm.R
#'
#' @import ggplot2
#'
#' @seealso [expectancy.HP()] and [expectancy.BLC()] for `HP` and `BLC` methods.
#'
#' [Heatmap.DLM()] and [Heatmap.list()] for `DLM` or `list` methods to drawing a Heatmap for the truncated life expectancy.
#'
#' @export
expectancy.DLM <- function(x, age = seq(0, max(fit$info$ages), by = 10),
                           graph = TRUE,
                           max_age = 110,
                           prob = 0.95,
                           ...){

  fit = x
  if(max(age) > max_age){
    stop("Invalid age interval. Check the max_age argument")
  }
  max_age = max_age+1

  #calculating qx and ci
  if(max_age > max(fit$info$ages)){
    pred <- predict( fit, h = (max_age - max(fit$info$ages)), prob = prob )
    mu <- c(fitted(fit)$qx_fitted,pred$qx_fitted)
    ic <- rbind(qx_ci(fit, prob = prob)[,-1],data.frame(qi = pred$qx_inf, qs = pred$qx_sup))
  }else{
    mu <- fitted(fit)$qx_fitted
    ic <- qx_ci(fit, prob=prob)
  }

  exp_total <- rep(NA, max_age)

  #cumprod for life expectancy (px)
  for (i in 1:max_age) {
    exp_total[i] <- sum(cumprod(1-mu[i:max_age]))
  }
  exp_total <- round(exp_total,2)


  exp_inf <- rep(NA,max_age); exp_sup <- rep(NA,max_age)

  ### upper CI:
  for (i in 1:max_age) {
    exp_sup[i] <- sum(cumprod(1-ic$qi[i:max_age]))
  }
  exp_sup <- round(exp_sup,2)

  ### lower CI:
  for (i in 1:max_age) {
    exp_inf[i] <- sum(cumprod(1-ic$qs[i:max_age]))
  }
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
    return(list(tabua=tab[tab$Age %in% age,],
                plot=p))
  }else{
    return(tab[tab$Age %in% age,])
  }
}


#' @export
expectancy.ClosedDLM <- function(x, age = seq(0, max(fit$info$ages), by = 10),
                                 graph = TRUE, prob = 0.95, ...){

  fit = x
  ###sanity check
  if(max(age) > max(fit$info$ages)){
    stop("Invalid age interval. Check the ages modeled")
  }

  #last age modeled
  max_age <- max(fit$info$ages)

  ##calculating log(qx)
  mu <- apply(fit$qx, 2, median)
  exp_total <- rep(NA, max_age)

  #cumprod for life expectancy (px)
  for (i in 1:max_age) {
    exp_total[i] <- sum(cumprod(1-mu[i:max_age]))
  }
  exp_total <- round(exp_total,2)


  ##ci
  ic <- qx_ci(fit, prob=prob)
  exp_inf <- rep(NA,max_age); exp_sup <- rep(NA,max_age)

  ### upper CI:
  for (i in 1:max_age) {
    exp_sup[i] <- sum(cumprod(1-ic$qi[i:max_age]))
  }
  exp_sup <- round(exp_sup,2)

  ### lower CI:
  for (i in 1:max_age) {
    exp_inf[i] <- sum(cumprod(1-ic$qs[i:max_age]))
  }
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
