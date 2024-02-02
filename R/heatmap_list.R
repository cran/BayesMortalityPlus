#' @name Heatmap.list
#' @rdname Heatmap.list
#'
#' @title Heatmap for a set of life tables
#'
#' @description This function plots a heatmap for the life expectancy of the mortality graduations
#'  returned by hp(), dlm(), hp_close() or dlm_close() functions.
#'
#'
#' @param x List of objects of classes: `HP`, `DLM`, `ClosedHP`, or `ClosedDLM`.
#' @param x_lab Description of the object 'fit'.
#' @param age Vector with the ages to plot the heatmap.
#' @param max_age Positive number indicating the last age to be considered to compute the life expectancy (methods for matching the age interval will be considered if needed). This argument is only necessary with objects of the class `HP` or `DLM`.
#' @param color Vector with colours used in the heatmap.
#' @param ... Further arguments passed to or from other methods.
#'
#' @return A ggplot2 heatmap of the life expectancy.
#'
#' @examples
#' ## Importing mortality data from the USA available on the Human Mortality Database (HMD):
#' data(USA)
#'
#' # Example (HP): -------------------------------
#'
#' ## Selecting the data from 2010
#' USA2010 = USA[USA$Year == 2010,]
#'
#' ExF = USA2010$Ex.Female[1:91]
#' DxF = USA2010$Dx.Female[1:91]
#' x <- 0:90
#'
#' fitF <- hp(x, ExF, DxF, model = "lognormal", M = 1000, bn = 0, thin = 10)
#'
#' ExM = USA2010$Ex.Male[1:91]
#' DxM = USA2010$Dx.Male[1:91]
#'
#' fitM <- hp(x, ExM, DxM, model = "lognormal", M = 1000, bn = 0, thin = 10)
#'
#' fits <- list(fitF = fitF, fitM = fitM)
#'
#' Heatmap(fits, x_lab = c("Female 2010 USA","Male 2010 USA"),
#'         age = 15:85)
#'
#'
#' @include expectancy_dlm.R
#' @include expectancy_hp.R
#'
#' @seealso [Heatmap.HP()], [Heatmap.DLM()] and [Heatmap.BLC()] for drawing single Heatmaps.
#'
#' @import ggplot2
#' @export
#'
Heatmap.list <- function(x, x_lab = NULL, age = NULL, max_age = NULL,
                         color = c("red","white","blue"), ...){

   fits = x
   if(is.null(x_lab)){
     x_lab <- rep(NA_character_,length(fits))
     for(i in 1:length(x_lab)){x_lab[i] <- paste("Fit",as.character(i))}
    }
   #sanity check:
   if(inherits(fits, "list")){
     if(length(fits) != length(x_lab)){stop("Number of fitted models is different of the x_lab's length.")}
   }
   if(length(color) != 3){stop("The argument color must be a 3 length vector.")}


   #checking the model
  if(is.null(age)){
    check = unlist(lapply(fits, class)) %in% c("DLM", "ClosedDLM")
    if(all(check)){
      ages = rep(NA_real_, length(check))
      for(i in 1:length(check)){ages[i] = length(fits[[i]]$info$ages)}
      age = fits[[which.min(ages)]]$info$ages
    }else if(!any(check)){
      age = 0:90
    }else{
      dlm_id <- seq(1:length(check))[check]
      ages = rep(NA_real_, length(dlm_id))
      for(i in dlm_id){ages[i] = length(fits[[i]]$info$ages)}
      age = fits[[which.min(ages)]]$info$ages
    }
  }

   #calculating life expectancy
   if( any(unlist(lapply(fits, class)) %in% c("ClosedHP","ClosedDLM")) ){
     if(!is.null(max_age)){
       warning("max_age argument is available for HP and DLM objects only, the argument will be dropped")
     }
     lista_exp <- lapply(fits, expectancy, graph = FALSE, age = age)
     exps = NULL
     for(i in 1:length(lista_exp)){exps <- rbind(exps,lista_exp[[i]])}
   }else{
     if(is.null(max_age)){max_age = 110}
     lista_exp <- lapply(fits, expectancy, graph = FALSE, age = age, max_age = max_age)
     exps = NULL
     for(i in 1:length(lista_exp)){exps <- rbind(exps,lista_exp[[i]])}
   }
   #creating heatmap dataframe
   exp <- exps$expectancy
   ano <- c()
   for(i in 1:length(x_lab)){ano <- c((rep(x_lab[i],length(age))),ano)}
   idade <- exps$age
   df <- data.frame(
     "age" = idade,
     "year" = rev(as.character(ano)),
     "exp" = exp)

   #plot
   midp <- mean(exp)
   p <- ggplot(df) + theme_light() +
     geom_raster(aes(x = year, y = age, fill = exp),interpolate = FALSE) +
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
