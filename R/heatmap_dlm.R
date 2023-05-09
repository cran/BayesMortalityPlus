#' @name Heatmap.DLM
#' @rdname Heatmap.DLM
#'
#' @title DLM: Heatmap for the life expectancy.
#'
#' @description This function plots a heatmap for the life expectancy of the fitted DLMs.
#'
#'
#' @param x Object or a list of objects of the class `DLM` or `ClosedDLM` returned by dlm() or dlm_close() functions.
#' @param x_lab Description of the object 'fit'.
#' @param age Vector with the ages to plot the heatmap.
#' @param max_age Positive number indicating the last age to be considered to compute the life expectancy (prediction will be considered to match the age interval if needed). This argument is only necessary with objects of the class `DLM`.
#' @param color Vector of colours used in the heatmap.
#' @param ... Further arguments passed to or from other methods.
#'
#' @return A ggplot2 heatmap of the life expectancy.
#'
#' @examples
#' ## Importing mortality data from the USA available on the Human Mortality Database (HMD):
#' data(USA)
#'
#' # Example 1: -------------------------------
#'
#' USA2010 = USA[USA$Year == 2010,]
#'
#' ExF = USA2010$Ex.Female[1:91]
#' DxF = USA2010$Dx.Female[1:91]
#'
#' qx_tF <- DxF/ExF
#' qx_tF <- 1 - exp(-qx_tF)
#' yF <- log(qx_tF)
#'
#' fitF <- dlm(yF, M = 100, bn = 20, thin = 1)
#'
#' Heatmap(fitF, x_lab = "Female expec. 2010 USA", max_age = 90)
#'
#'
#' @include expectancy_dlm.R
#'
#' @seealso [Heatmap.BLC()] and [Heatmap.HP()] for `BLC` or `HP` methods.
#'
#' [Heatmap.list()] to the `list` method, adding multiple objects in one single Heatmap.
#'
#' @import ggplot2
#' @export
#'
Heatmap.DLM <- function(x, x_lab, age = NULL, max_age = 110,
                        color = c("red","white","blue"), ...){
  fits = x
  if(is.null(age)){age = fits$info$ages}
  #sanity check:
  if(max(age) > max_age){
    stop("Invalid age interval. Check the max_age argument")
  }

  if(length(color) != 3){stop("The argument color must be a 3 length vector.")}


  #calculating the life expectancy
  exps <- expectancy(fits, graph = FALSE, age = age, max_age = max_age)

  #creating dataframe for the heatmap:
  exp <- exps$Expectancy
  ano <- c()
  for(i in 1:length(x_lab)){
    ano <- c((rep(x_lab[i],length(age))),ano)
  }
  idade <- exps$Age
  df <- data.frame(
    "age" = idade,
    "year" = rev(as.character(ano)),
    "exp" = exp)

  #plot
  midp <- mean(exp)

  p <- ggplot(df) + theme_light() +
    geom_raster(aes(x = year, y = age, fill = exp),interpolate = TRUE) +
    labs(x="",
         y="Age",
         title = "Life expectancy") +
    scale_fill_gradient2(name = "Expectancy (years)",
                         low = color[1],
                         mid = color[2],
                         high = color[3],
                         midpoint = midp)

  return(p)
}


#' @export
Heatmap.ClosedDLM <- function(x, x_lab, age = NULL,
                              color = c("red","white","blue"), ...){
  fits = x
  #sanity check:
  if(is.null(age)){age = fits$info$ages}
  if(max(age) > max(fits$info$ages)){
    stop("Invalid age interval. Check the max_age argument")
  }

  if(length(color) != 3){stop("The argument color must be a 3 length vector.")}

  #calculating the life expectancy
  exps <- expectancy(fits, graph = FALSE, age = age)

  #creating dataframe for the heatmap:
  exp <- exps$Expectancy
  idade <- exps$Age
  ano <- c()
  for(i in 1:length(x_lab)){
    ano <- c((rep(x_lab[i],length(age))),ano)
  }
  df <- data.frame(
    "age" = idade,
    "year" = rev(as.character(ano)),
    "exp" = exp)

  #plot
  midp <- mean(exp)

  p <- ggplot(df) + theme_light() +
    geom_raster(aes(x = year, y = age, fill = exp),interpolate = TRUE) +
    labs(x="",
         y="Age",
         title = "Life expectancy") +
    scale_fill_gradient2(name = "Expectancy (years)",
                         low = color[1],
                         mid = color[2],
                         high = color[3],
                         midpoint = midp)

  return(p)
}

