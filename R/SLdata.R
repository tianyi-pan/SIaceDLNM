#' Prepare Super Learner quantities for an SIaceDLNM fit
#'
#' @param object An object of class \code{acmeDLNM_fit}.
#' @param dat.c Required counterfactual data, with the same reference records
#'   and full exposure histories, on the predictor scale used by \code{predict()}.
#' @param verbose Whether to print progress, default \code{FALSE}.
#' @param AIC.only Prepare singleton AIC quantities instead of cross-validation
#'   predicted means. Used by \code{mDLNM::AICdata()}.
#' @param nthreads Number of threads for deletion calculations, default \code{1}.
#'   Uses OpenMP for parallel calculation.
#' @param ... Unused; additional arguments are rejected.
#' @exportS3Method mDLNM::SLdata
SLdata.acmeDLNM_fit <- function(object, dat.c, verbose = FALSE,
                                AIC.only = FALSE, nthreads = 1, ...) {
  if(missing(dat.c) || is.null(dat.c)) stop("dat.c is required.")
  if(length(list(...))) stop("Unused arguments.")
  if(identical(object$family, "gaussian")) stop("SLdata() requires a negative binomial fit.")

  dat <- object$data
  point <- object$point
  coef.samples <- NULL
  if(AIC.only) {
    samples <- object$CI.sample
    if(is.null(samples$alpha_f_sample) || nrow(samples$alpha_f_sample) < 2) {
      stop("AICdata() requires at least two stored coefficient samples in CI.sample.")
    }
    index.coef <- t(solve(dat$Bindex, t(samples$index_sample)))
    index.coef <- index.coef[, -1, drop = FALSE] / index.coef[, 1]
    coef.samples <- cbind(samples$alpha_f_sample, samples$phi_sample, index.coef,
                          samples$betaR_sample, samples$betaF_sample)
  }
  mod.build <- acmeDLNMbuild(dat$y, dat$B_inner, dat$knots_f,
                            dat$SwI_large, dat$SfI, dat$Dw,
                            dat$Xrand, dat$Xfix, dat$Zf %*% dat$Ufpen,
                            dat$offset$Xoffset, dat$r, dat$K, dat$a, dat$Bindex,
                            point$alpha_f, point$phi, point$delta,
                            point$log_theta, point$log_smoothing_f, point$log_smoothing_w,
                            point$betaR, point$betaF, point$log_smoothing)
  ptr.cf <- stats::predict(object, dat = dat.c, .build.only = TRUE)$address.eigen
  out <- SLacmeDLNM(mod.build$address.eigen, ptr.cf, verbose, AIC.only,
                   nthreads, coef.samples)
  theta <- exp(point$log_theta)

  if(AIC.only) {
    out$aic <- as.vector(out$aic)
    names(out$aic) <- colnames(out$aic.loo) <- c("AIC_conventional", "AIC_cAIC")
    return(list(y = dat$y, row.id = rownames(object$modeldata), theta = theta,
                aic = out$aic, aic.loo = out$aic.loo,
                delta = out$delta, delta.loo = as.vector(out$delta.loo),
                variance = stats::var(as.vector(out$delta.sample)),
                estimand = "ratio"))
  }

  return(list(y = dat$y, row.id = rownames(object$modeldata), theta = theta,
              mu.loo = as.vector(out$mu.loo), mu.lto = out$mu.lto,
              delta = out$delta, delta.loo = as.vector(out$delta.loo),
              estimand = "ratio"))
}
