#' Title Rate Ratio for acmeDLNM_fit objects
#'
#' @param object object of class \code{acmeDLNM_fit}.
#' @param verbose whether to print messages during the process, default \code{FALSE}.
#' @details Rate ratios are only defined here for negative binomial fits with
#'   a log link. This method does not support Gaussian fits with an identity link.
#'
#' @return
#' @importFrom mgcv s
#' @exportS3Method mDLNM::RateRatio
RateRatio.acmeDLNM_fit <- function(object, x0, x1, verbose = FALSE, ...) {
  if(identical(object$family, "gaussian")) {
    stop("RateRatio() is not available for Gaussian fits with an identity link.")
  }
  ## point estimate
  pc <- object$data$pc
  alpha_w <- object$point$alpha_w
  alpha_f <- object$point$alpha_f
  kw <- object$data$kw
  kE <- object$data$kE

  Ufpen <- object$data$Ufpen
  Uwpen <- object$data$Uwpen
  Zwnew <- object$data$Zw %*% Uwpen
  maxL <- object$data$maxL


  ## point estimate
  wl.fit <- function(lnew) c(1, Bsplinevec2Con(lnew, object$data$knots_w, 4, Zwnew)) %*% alpha_w
  wl.fit <- Vectorize(wl.fit)

  wl.discrete <- sapply(0:maxL, function(ll) {
    integrate(wl.fit, lower = ll, upper = ll + 1)$value
  })
  E0 <- sum(wl.discrete * Reduce("+", mapply("*", x0, object$point$index, SIMPLIFY = FALSE)))
  E1 <- sum(wl.discrete * Reduce("+", mapply("*", x1, object$point$index, SIMPLIFY = FALSE)))

  fE.fit <- function(Enew) as.numeric(mgcv::PredictMat(object$smooth$fE, data = data.frame(E = Enew)) %*% Ufpen %*% alpha_f)
  eta.E.0 <- fE.fit(E0)
  eta.E.1 <- fE.fit(E1)

  RR.est <- data.frame(eta.E.0 = eta.E.0,
                       eta.E.1 = eta.E.1,
                       RR = exp(eta.E.1 - eta.E.0))

  ## CI
  RR.sample <- lapply(1:nrow(object$CI.sample$alpha_w_sample), function(i) {
    if(verbose) {
      if(i %% 100 == 0) {
        cat("Sampling iteration:", i, " / ", nrow(object$CI.sample$alpha_w_sample), "\n")
      }
    }
    wl.fit <- function(lnew) c(1, Bsplinevec2Con(lnew, object$data$knots_w, 4, Zwnew)) %*% object$CI.sample$alpha_w_sample[i,]
    wl.fit <- Vectorize(wl.fit)

    wl.discrete <- sapply(0:maxL, function(ll) {
      integrate(wl.fit, lower = ll, upper = ll + 1)$value
    })
    E0 <- sum(wl.discrete * Reduce("+", mapply("*", x0, object$CI.sample$index[i,], SIMPLIFY = FALSE)))
    E1 <- sum(wl.discrete * Reduce("+", mapply("*", x1, object$CI.sample$index[i,], SIMPLIFY = FALSE)))

    fE.fit <- function(Enew) as.numeric(mgcv::PredictMat(object$smooth$fE, data = data.frame(E = Enew)) %*% Ufpen %*% object$CI.sample$alpha_f_sample[i,])
    eta.E.0 <- fE.fit(E0)
    eta.E.1 <- fE.fit(E1)
    return(data.frame(eta.E.0 = eta.E.0,
                      eta.E.1 = eta.E.1,
                      RR = exp(eta.E.1 - eta.E.0)))
  })


  RR.sample <- data.table::rbindlist(RR.sample)
  return(list(RR.est = RR.est,
              RR.sample = RR.sample)
         )
}
