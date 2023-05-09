#' @name fitted.BLC
#' @rdname fitted.BLC
#'
#' @title Fitted values for BLC fitted models
#'
#' @description Computes the fitted values associated to each age and year based on
#' the resulting chains from a fitted BLC model. In addition, this function also
#' evaluates the values of lower and upper limits of the credible interval.
#'
#'
#' @param object A `BLC` object, result of a call to blc() function.
#' @param cred A real number that indicates the probability of the credible interval.
#' @param ... Other arguments.
#'
#' @return A list with the matrices of fitted values and lower and upper limits of the credible interval for each age and year.
#'
#' @examples
#' ## Importing log-mortality data from Portugal:
#' data(PT)
#' Y <- PT
#'
#' ## Fitting the model
#' fit = blc(Y = Y, numit = 100, warmup = 20)
#'
#' ## Log-mortalities estimates for each age and year in model fitted
#' fitted(fit, cred = 0.95)
#'
#' @seealso [fitted.HP()] and [fitted.DLM()] for `HP` or `DLM` methods.
#'
#' @export
fitted.BLC <- function(object, cred = 0.95, ...) {
  obj = object
	N <- obj$numit - obj$warmup
	L <- nrow(obj$kappa)
	q <- nrow(obj$alpha)
	fits <- array(dim = c(q, L, N))

	for (i in (obj$warmup+1):obj$numit) {
		fits[ , ,i - obj$warmup] <- obj$alpha[ ,i] + obj$beta[ ,i,drop=F] %*% obj$kappa[ ,i]
	}

	alpha <- 1 - cred
	mean <- apply(fits, c(1,2), mean)
	upper <- apply(fits, c(1,2), quantile, 1 - alpha/2)
	lower <- apply(fits, c(1,2), quantile, alpha/2)

	colnames(mean) <- colnames(obj$Y)
	colnames(upper) <- colnames(obj$Y)
	colnames(lower) <- colnames(obj$Y)

	row.names(mean) <- row.names(obj$Y)
	row.names(upper) <- row.names(obj$Y)
	row.names(lower) <- row.names(obj$Y)

	list(mean = exp(mean), upper = exp(upper), lower = exp(lower))
}



#'
#' @export
fitted.PredBLC <- function(object, cred = 0.95, ...) {
  obj = object
  fits <- obj$y


  alpha <- 1 - cred
  mean <- apply(fits, c(3,2), mean)
  upper <- apply(fits, c(3,2), quantile, 1 - alpha/2)
  lower <- apply(fits, c(3,2), quantile, alpha/2)

  colnames(mean) <- colnames(obj$y)
  colnames(upper) <- colnames(obj$y)
  colnames(lower) <- colnames(obj$y)

  row.names(mean) <- row.names(obj$y)
  row.names(upper) <- row.names(obj$y)
  row.names(lower) <- row.names(obj$y)

  list(mean = exp(mean), upper = exp(upper), lower = exp(lower))
}
