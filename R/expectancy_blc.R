#' @name expectancy.BLC
#' @rdname expectancy.BLC
#'
#' @title Life Expectancy for fitted models and forecast
#'
#' @description Computes the fitted life expectancy for a specific age for each year in fit or prediction. It also calculates the limits of credible intervals.
#'
#' @param x A `BLC` or `PredBLC` object.
#' @param at A number that determines at which age the expectancy life is calculated based on the ages used in fit or prediction. For instance, at = 0 is related to the first age used in fitted model.
#' @param cred A number that specifies the probability of the credible interval. Default is '0.95'.
#' @param ... Further arguments passed to or from other methods.
#'
#' @return A list that contains three vectors with the fitted values of life expectancy and the lower and upper limits of the credible intervals for each year used in fitted model or for the prediction.
#'
#' @examples
#' ## Importing log-mortality data from Portugal:
#' data(PT)
#' Y <- PT
#'
#' ## Fitting the model
#' fit = blc(Y = Y, numit = 100, warmup = 20)
#'
#' ## Life expectancy for the years used in model fitted
#' expectancy(fit)
#'
#' ## Life expectancy for the tenth and thirtieth age in the years used in
#' ## model fitted (27 and 47 y.o.)
#' expectancy(fit, at = c(10,30))
#'
#' @seealso [expectancy.HP()] and [expectancy.DLM()] for `HP` and `DLM` methods.
#'
#' [Heatmap.BLC()] for `BLC` method to drawing a Heatmap for the truncated life expectancy.
#'
#' @export
expectancy.BLC <- function(x, at = NULL, cred = 0.95, ...) {
  obj <- x
	objClass <- class(obj)
	supportedClasses <- c("BLC", "ARBLC")

	if (!any(objClass %in% supportedClasses)) {
		stop("Invalid object type")
	}

	# N <- obj$numit - obj$warmup
	L <- nrow(obj$kappa)
	q <- nrow(obj$alpha)

	if (!is.null(at)){
	  if (mode(at) != "numeric"){stop("Expected `at` to be numeric")}
	  if (min(at) == 0){stop("Expected `at` to be greater than 0")}
	}else{
	  at <- 1:q
	}

	## Fitted values for observed series
	# fits <- array(dim = c(q, L, N))
	#
	# for (i in (obj$warmup+1):obj$numit) {
	# 	fits[ , ,i - obj$warmup] <- obj$alpha[ ,i] + obj$beta[ ,i,drop=F] %*% obj$kappa[ ,i]
	# }

	fits <- fitted(obj, cred = cred)

	exp_total <- matrix(NA, nrow = q, ncol = L)

	#cumprod for life expectancy (px)
	for (i in 1:(q-1)){
	  exp_total[i,] <- apply(fits$mean[i:q,], 2, function(x) sum(cumprod(1-x)))
	}
	exp_total[q,] <- 1 - (fits$mean[q,])
	exp_total <- round(exp_total,3)

	##ci
	exp_inf <- matrix(NA, nrow = q, ncol = L); exp_sup <- matrix(NA, nrow = q, ncol = L)

	### upper CI:
	for (i in 1:(q-1)){
	  exp_sup[i,] <- apply(fits$lower[i:q,], 2, function(x) sum(cumprod(1-x)))
	}
	exp_sup[q,] <- 1 - (fits$lower[q,])
	exp_sup <- round(exp_sup,3)

	### lower CI:
	for (i in 1:(q-1)){
	  exp_inf[i,] <- apply(fits$upper[i:q,], 2, function(x) sum(cumprod(1-x)))
	}
	exp_inf[q,] <- 1 - (fits$upper[q,])
	exp_inf <- round(exp_inf,3)

   ### IMPLEMENTAR SAIDA MAIS PARECIDA COM HP E DLM, COM DATA FRAMES COM IDADES SELECIONADAS
	## PARA TODOS OS ANOS, LISTA COM $MEAN $LOWER $UPPER
	colnames(exp_total) <- colnames(obj$Y)
	colnames(exp_sup) <- colnames(obj$Y)
	colnames(exp_inf) <- colnames(obj$Y)

	row.names(exp_total) <- row.names(obj$Y)
	row.names(exp_sup) <- row.names(obj$Y)
	row.names(exp_inf) <- row.names(obj$Y)

	ret <- list(expectancy = exp_total[at,], upper = exp_sup[at,], lower = exp_inf[at,])

	ret
}


#'
#' @export
expectancy.PredBLC <- function(x, at = NULL, cred = 0.95, ...) {
  obj <- x

  if (!inherits(obj,"PredBLC")) stop("Invalid object")

  L <- obj$h
  q <- dim(obj$y)[3]

  if (!is.null(at)){
    if (mode(at) != "numeric"){stop("Expected `at` to be numeric")}
    if (min(at) == 0){stop("Expected `at` to be greater than 0")}
  }else{
    at <- 1:q
  }

  fits <- fitted(obj, cred = cred)

  exp_total <- matrix(NA, nrow = q, ncol = L)

  for (i in 1:(q-1)){
    exp_total[i,] <- apply(fits$mean[i:q,], 2, function(x) sum(cumprod(1-x)))
  }
  exp_total[q,] <- 1 - (fits$mean[q,])
  exp_total <- round(exp_total,3)

  # alpha <- 1 - cred
  # h <- dim(obj$y)[2]
  #
  # ages <- (at):(dim(obj$y)[3]-1)
  # ret <- list(expectancy = numeric(h), lower = numeric(h), upper = numeric(h))
  #
  # for (i in 1:h) {
  #   tmp <- apply(obj$y[ ,i, ages+1], 1, function(x) sum(cumprod(1-exp(x))))
  #   #somatorio da multiplicacao acumulada das exp de vida truncadas
  #   ret$expectancy[i] <- mean(tmp)
  #   ret$lower[i] <- quantile(tmp, alpha/2)
  #   ret$upper[i] <- quantile(tmp, 1 - alpha/2)
  # }

  ##ci
  exp_inf <- matrix(NA, nrow = q, ncol = L); exp_sup <- matrix(NA, nrow = q, ncol = L)

  ### upper CI:
  for (i in 1:(q-1)){
    exp_sup[i,] <- apply(fits$lower[i:q,], 2, function(x) sum(cumprod(1-x)))
  }
  exp_sup[q,] <- 1 - (fits$lower[q,])
  exp_sup <- round(exp_sup,3)

  ### lower CI:
  for (i in 1:(q-1)){
    exp_inf[i,] <- apply(fits$upper[i:q,], 2, function(x) sum(cumprod(1-x)))
  }
  exp_inf[q,] <- 1 - (fits$upper[q,])
  exp_inf <- round(exp_inf,3)

  ret <- list(expectancy = exp_total[at,], upper = exp_sup[at,], lower = exp_inf[at,])

  ret
}
