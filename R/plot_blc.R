#' @name plot.BLC
#' @rdname plot.BLC
#'
#' @title BLC: Plot the fitted values
#'
#' @description This function plots the fitted log mortality values as well as the parameters values and credible intervals of the BLC fitted models.
#'
#'
#' @param x A `BLC` object, result of a call to blc() function.
#' @param parameter A character determines the parameter that will be plotted. Default is "all" which means that all three parameters "alpha", "beta" and "kappa" will be plotted. It can also be "alpha", "beta", "kappa" or "fitted". The last one provides a plot with all the fitted tables.
#' @param prob A numeric value that indicates the probability for the credible interval. Default is '0.9'.
#' @param age A numeric vector that represents the ages used in the fitted BLC model. Default is 'NULL'.
#' @param ... Other arguments.
#'
#' @return A plot with the fitted log mortality or fitted values and credible intervals of the parameters.
#'
#' @examples
#' ## Importing log-mortality data from Portugal:
#' data(PT)
#' Y <- PT
#'
#' ## Fitting the model
#' fit = blc(Y = Y, M = 100, bn = 20)
#'
#' ## Parameters' plot
#' plot(fit, parameter = "all")
#' \donttest{plot(fit, parameter = "beta", prob = 0.95)
#' plot(fit, parameter = "alpha", age = 18:80)
#' plot(fit, parameter = "kappa")
#'
#' ## Fitted mortality graduation
#' plot(fit, parameter = "fitted", age = 18:80)
#' }
#'
#' @seealso [plot.HP()] and [plot.DLM()] for `HP` or `DLM` methods.
#'
#' @import ggplot2
#' @import scales
#'
#' @export
plot.BLC <- function(x, parameter = "all", prob = 0.9,
                     age = NULL, ...) {
  obj = x
  sig = 1 - prob
  if(!(is.null(age))){
    if(!(is.integer(age))){
      stop("Object age should be an integer vector")
    }
  }else{
    age = 1:length(obj$beta[,1])
  }

  # Drop warmup from chains
  chain.idx <- obj$bn:obj$M
  alpha <- obj$alpha[ ,chain.idx]
  beta <- obj$beta[ ,chain.idx]
  kappa <- obj$kappa[ ,chain.idx]

  alpha.est <- apply(alpha, 1, mean)
  alpha.inf <- apply(alpha, 1, quantile, sig/2)
  alpha.sup <- apply(alpha, 1, quantile, 1 - sig/2)
  alpha.lim <- range(c(alpha.inf, alpha.sup))

  beta.est <- apply(beta, 1, mean)
  beta.inf <- apply(beta, 1, quantile, sig/2)
  beta.sup <- apply(beta, 1, quantile, 1 - sig/2)
  beta.lim <- range(c(beta.inf, beta.sup))

  kappa.est <- apply(kappa, 1, mean)
  kappa.inf <- apply(kappa, 1, quantile, sig/2)
  kappa.sup <- apply(kappa, 1, quantile, 1 - sig/2)
  kappa.lim <- range(c(kappa.inf, kappa.sup))

  N <- length(kappa.est)

  # Plot types
  if(parameter == "all"){

    df.beta = data.frame(x = age, fitted = beta.est, lim.inf = beta.inf, lim.sup = beta.sup, param = "beta")
    df.alpha = data.frame(x = age, fitted = alpha.est, lim.inf = alpha.inf, lim.sup = alpha.sup, param = "alpha")
    df.kappa = data.frame(x = 1:N, fitted = kappa.est, lim.inf = kappa.inf, lim.sup = kappa.sup, param = "kappa")
    df = rbind(df.beta, df.alpha, df.kappa)
    df$param = factor(df$param, labels = c("alpha[x]", "beta[x]", "kappa[t]"))

    ggplot(data = df) +
      geom_ribbon(aes(x = x, ymin = lim.inf, ymax = lim.sup), alpha = 0.5, fill = "blue") +
      geom_line(aes(x = x, y = fitted), col = "blue") +
      xlab("") + ylab("") + theme_bw() +
      theme(axis.title.x = ggplot2::element_text(color = 'black', size = 13),
            axis.title.y = ggplot2::element_text(color = 'black', size = 13)) +
      facet_wrap(~param, scales = "free", nrow = 2,
                 labeller = label_parsed) +
      geom_hline(data = data.frame(yint = 0, param = "beta[x]"), aes(yintercept = yint), col = "red", lty = 2)

  }else if(parameter == "beta"){

    ggplot(data = data.frame(x = age, fitted = beta.est, lim.inf = beta.inf, lim.sup = beta.sup)) +
      geom_hline(yintercept = 0, col = "red", lty = 2) +
      geom_ribbon(aes(x = x, ymin = lim.inf, ymax = lim.sup), alpha = 0.5, fill = "blue") +
      geom_line(aes(x = x, y = fitted), col = "blue") +
      xlab("x") + ylab("qx") + theme_bw() +
      theme(axis.title.x = ggplot2::element_text(color = 'black', size = 13),
            axis.title.y = ggplot2::element_text(color = 'black', size = 13))

  }else if(parameter == "alpha"){

    ggplot(data = data.frame(x = age, fitted = alpha.est, lim.inf = alpha.inf, lim.sup = alpha.sup)) +
      geom_ribbon(aes(x = x, ymin = lim.inf, ymax = lim.sup), alpha = 0.5, fill = "blue") +
      geom_line(aes(x = x, y = fitted), col = "blue") +
      xlab("x") + ylab(expression(alpha[x])) + theme_bw() +
      theme(axis.title.x = ggplot2::element_text(color = 'black', size = 13),
            axis.title.y = ggplot2::element_text(color = 'black', size = 13))

  }else if(parameter == "kappa"){

    N = length(kappa.est)
    ggplot(data = data.frame(x = 1:N, fitted = kappa.est, lim.inf = kappa.inf, lim.sup = kappa.sup)) +
      geom_ribbon(aes(x = x, ymin = lim.inf, ymax = lim.sup), alpha = 0.5, fill = "blue") +
      geom_line(aes(x = x, y = fitted), col = "blue") +
      xlab("t") + ylab(expression(kappa[t])) + theme_bw() +
      theme(axis.title.x = ggplot2::element_text(color = 'black', size = 13),
            axis.title.y = ggplot2::element_text(color = 'black', size = 13))

  }else if(parameter == "fitted"){

    N = length(age); t = length(kappa.est)
    tables = matrix(NA_real_, nrow = N, ncol = t)
    aux.table = matrix(NA_real_, nrow = N, ncol = ncol(alpha))
    for(t in 1:t){
      for(i in 1:ncol(alpha)){  aux.table[,i] = alpha[,i] + beta[,i]*kappa[t,i] }
      tables[,t] = rowMeans(aux.table)
    }
    df.tables = data.frame(tables, Age = age) %>% gather(key = "Year", value = "log.qx", -Age)
    ggplot(df.tables) +
      scale_y_continuous(trans = "log10", breaks = 10^-seq(0,5), limits = 10^-c(5,0), labels = scales::comma) +
      scale_x_continuous(breaks = seq(0, 100, by = 10), limits = c(NA,NA)) +
      geom_line(aes(x = Age, y = exp(log.qx), col = Year), show.legend = FALSE) +
      xlab("x") + ylab("qx") + theme_bw() +
      theme(axis.title.x = ggplot2::element_text(color = 'black', size = 13),
            axis.title.y = ggplot2::element_text(color = 'black', size = 13))

  }
}

