#' @name Heatmap.BLC
#' @rdname Heatmap.BLC
#'
#' @title BLC: Heatmap for the life expectancy
#'
#' @description Draws a Heat Map based on the life expectancy of a fitted BLC or PredBLC model.
#'
#'
#' @param x A `BLC` or `PredBLC` object, result of a call to blc() function or forecast via predict() function.
#' @param x_lab Description of the modelled object.
#' @param age Vector with the ages to plot the heatmap.
#' @param color Vector of colours used in the heatmap.
#' @param ... Further arguments passed to or from other methods.
#'
#' @return A ggplot2 heatmap of the life expectancy.
#'
#' @examples
#' ## Importing log-mortality data from Portugal:
#' data(PT)
#' Y <- PT
#'
#' ## Fitting the model
#'
#' fit = blc(Y = Y, M = 100, bn = 20)
#'
#' ## Heatmap:
#' Heatmap(fit, x_lab = 2000:2015, age = 18:80)
#'
#' @seealso [Heatmap.HP()] and [Heatmap.DLM()] for `HP` or `DLM` methods.
#'
#' @import ggplot2
#' @export
Heatmap.BLC <- function(x, x_lab = NULL, age = NULL, color = c("red","white","blue"), ...){
  obj = x
  objClass <- class(obj)
  supportedClasses <- c("BLC")

  if(!any(objClass %in% supportedClasses)){stop("Invalid object type")}

  L <- ncol(obj$Y)
  q <- nrow(obj$Y)

  exps <- expectancy(obj)$expectancy

  if(!is.null(x_lab)){
    if(length(x_lab) != L){stop("Argument `x_lab` has to be the same length of the modelled years")}
    tag <- x_lab
  }else{tag <- 1:L}
  if(!is.null(age)){
    if(length(age) != q){stop("Argument `age` has to be the same length of the modelled age interval")}
    rows <- age
  }else{rows <- 1:q}

  exps2 <- as.data.frame(matrix(NA_real_, nrow = q*L, ncol = 3))
  colnames(exps2) <- c("Exp","Year","Age")
  for(j in 1:L){exps2[(q*j-(q-1)):(q*j),] <- data.frame(exps[,j],rep(tag[j]),as.numeric(rows))}
  midp <- mean(exps)

  p <- ggplot(exps2) + theme_light() +
    ggplot2::theme(plot.title = ggplot2::element_text(lineheight = 1.2),
                   axis.title.x = ggplot2::element_text(color = 'black', size = 12),
                   axis.title.y = ggplot2::element_text(color = 'black', size = 12),
                   axis.text = ggplot2::element_text(color = 'black', size = 12),
                   legend.text = ggplot2::element_text(size = 12),
                   legend.position = "bottom") +
  		geom_tile(aes(x = reorder(as.character(Year), sort(as.numeric(Year))), y = Age, fill=Exp)) +
  		labs(x="Years",
  		      y="Age") +
  	  scale_fill_gradient2(name = "Expectancy (years)",
  	                       low = color[1],
  	                       mid = color[2],
  	                       high = color[3],
  	                       midpoint = midp)

  p
}

#'
#' @export
Heatmap.PredBLC <- function(x, x_lab = NULL, age = NULL, color = c("red","white","blue"), ...){
  obj = x
  objClass <- class(obj)
  supportedClasses <- c("PredBLC")

  if(!any(objClass %in% supportedClasses)){stop("Invalid object type")}

  L <- obj$h
  q <- dim(obj$y)[3]

  exps <- expectancy(obj)$expectancy

  if(!is.null(x_lab)){
    if(length(x_lab) != L){stop("Argument `x_lab` has to be the same length of the years forecasted")}
    tag <- x_lab
    }else{tag <- 1:L}
  if(!is.null(age)){
    if(length(age) != q){stop("Argument `age` has to be the same length of the modelled age interval")}
    rows <- age
    }else{rows <- 1:q}

  exps2 <- as.data.frame(matrix(NA_real_, nrow = q*L, ncol = 3))
  colnames(exps2) <- c("Exp","Year","Age")
  for(j in 1:L){exps2[(q*j-(q-1)):(q*j),] <- data.frame(exps[,j], rep(tag[j]),as.numeric(rows))}
  midp <- mean(exps)

  p <- ggplot(exps2) + theme_light() +
    geom_tile(aes(x = reorder(as.character(Year), sort(Year)), y = Age, fill=Exp)) +
    labs(x="Years",
         y="Age",
         title = "Life expectancy") +
    scale_fill_gradient2(name = "Expectancy (years)",
                         low = color[1],
                         mid = color[2],
                         high = color[3],
                         midpoint = midp)

  p
}
