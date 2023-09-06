#' @name Heatmap.HP
#' @rdname Heatmap.HP
#'
#' @title HP: Heatmap for the life expectancy
#'
#' @description This function plots a heatmap for the life expectancy of the fitted HP models.
#'
#'
#' @param x Object or a list of objects of the class `HP` or `ClosedHP` returned by hp() or close_hp() functions.
#' @param x_lab Description of the object 'fit'.
#' @param age Vector with the ages to plot the heatmap.
#' @param max_age Positive number indicating the last age to be considered to compute the life expectancy (extrapolation will be considered to match the age interval if needed). This argument is only necessary with objects of the class `HP`.
#' @param color Vector of colours used in the heatmap.
#' @param ... Further arguments passed to or from other methods.
#'
#' @return A ggplot2 heatmap of the life expectancy.
#'
#' @examples
#' ## Importing mortality data from the USA available on the Human Mortality Database (HMD):
#' data(USA)
#'
#' # Example: -------------------------------
#'
#' USA2010 = USA[USA$Year == 2010,]
#'
#' ExF = USA2010$Ex.Female[1:91]
#' DxF = USA2010$Dx.Female[1:91]
#' x <- 0:90
#'
#' fitF <- hp(x, ExF, DxF, model = "lognormal", M = 1000, bn = 0, thin = 10)
#'
#' Heatmap(fitF, x_lab = "Female expec. 2010 USA")
#'
#' @seealso [Heatmap.BLC()] and [Heatmap.DLM()] for `BLC` or `DLM` methods.
#'
#' [Heatmap.list()] to the `list` method, adding multiple objects in one single Heatmap.
#'
#' @include expectancy_hp.R
#'
#' @import ggplot2
#' @export
#'
Heatmap.HP <- function(x, x_lab = NULL, age = 0:90, max_age = 110,
                       color = c("red","white","blue"), ...){
  fits = x
  if(is.null(x_lab)){x_lab <- "Fitted model"}
  #sanity check:
  if(max(age) > max_age){stop("Invalid age interval. Check the max_age argument")}
  if(length(color) != 3){stop("The argument color must be a 3 length vector.")}

  #calculating life expectancy
  exps <- expectancy(fits, graph = FALSE, age = age, max_age = max_age)


  #creating dataframe for the heatmap:
  exp <- exps$Expectancy
  ano <- c()
  for(i in 1:length(x_lab)){ano <- c((rep(x_lab[i],length(age))),ano)}
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
#'
Heatmap.ClosedHP <- function(x, x_lab = NULL, age = 0:90,
                             color = c("red","white","blue"), ...){
  fits = x
  if(is.null(x_lab)){x_lab <- "Fitted model"}

  #sanity check:
  if(length(color) != 3){stop("The argument color must be a 3 length vector.")}

  #calculating life expectancy
  exps <- expectancy(fits, graph = FALSE, age = age)

  #creating dataframe for the heatmap:
  exp <- exps$Expectancy
  ano <- c()
  for(i in 1:length(x_lab)){ano <- c((rep(x_lab[i],length(age))),ano)}
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
