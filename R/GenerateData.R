#' Title Generate data
#'
#' @param fEtype
#' @param wltype
#' @param Nt
#' @param kx.per500
#' @param theta
#' @param maxL
#' @param kw
#'
#' @return
#' @export
#' @importFrom mgcv s
#'
#' @examples
GenerateData <- function(fEtype, wltype, Nt = 1000, kx.per500 = 100,
                         interpolate = TRUE,
                         theta = 10,
                         M = 3,
                         gamma = c(1,0.5,0.5),
                         maxL = 14,
                         exposure = NULL,
                         std = TRUE,
                         verbose = FALSE,
                         c1 = NULL,
                         c2 = NULL,
                         wl.set = NULL,
                         fE.set = NULL,
                         other) {

  ## load waterloo PM2.5 data (PM25.waterloo)
  if(is.null(exposure)) data("AirPollutionWaterloo")
  maxLreal <- maxL+1

  t <- 1:(Nt+maxL) # time 1 to (Nt+40)
  t.sim <- t[-(1:maxL)]

  if((kx.per500 > 300) || (interpolate == TRUE)) {
    # cat("Interpolate the exposure process.")
    kx <- Nt+maxL + 4 + 2 # number of knots for X(t)
    interpolate <- TRUE
  } else{
    kx <- kx.per500 * ifelse(Nt < 500, 1, round(Nt/500)) # number of knots for X(t)
  }


  ## exposure process
  if(is.null(exposure)) {
    PM25 <- AirPollution.Waterloo$PM25[t]
    O3 <- AirPollution.Waterloo$O3[t]
    NO2 <- AirPollution.Waterloo$NO2[t]
  }

  PM25.orig <- PM25
  O3.orig <- O3
  NO2.orig <- NO2

  if(std) {
    PM25 <- (PM25 - mean(PM25)) / sd(PM25)
    O3 <- (O3 - mean(O3)) / sd(O3)
    NO2 <- (NO2 - mean(NO2)) / sd(NO2)
  }

  if(is.null(wl.set)) {
    switch (wltype,
            type1 = {
              wlc <- function(l) dnorm(l, mean = 3, sd = 3.5)^2
              wl_de <- sqrt(integrate(wlc, lower = 0, upper = maxLreal)$value)
              wl <- function(l) dnorm(l, mean = 3, sd = 3.5)/wl_de
            },
            type2 = {
              wlc <- function(l) (1/(1+exp(l-8)))^2
              wl_de <- sqrt(integrate(wlc, lower = 0, upper = maxLreal)$value)
              wl <- function(l) 1/(1+exp(l-8))/wl_de
            },
            type3 = {
              wlc <- function(l) (1/(1+exp(l-1)))^2
              wl_de <- sqrt(integrate(wlc, lower = 0, upper = maxLreal)$value)
              wl <- function(l) 1/(1+exp(l-1))/wl_de
            },
            type4 = {
              wlc <- function(l) dnorm(l,8,10)^2
              wl_de <- sqrt(integrate(wlc, lower = 0, upper = maxLreal)$value)
              wl <- function(l) dnorm(l,8,10)/wl_de
            }
            # symmetric = {
            #   wlc <- function(l) ( dnorm(l, mean = 7.5, sd = 3.2)-dnorm(0, mean = 7.5, sd = 3.2) )^2
            #   wl_de <- sqrt(integrate(wlc, lower = 0, upper = maxLreal)$value)
            #   wl <- function(l) ( dnorm(l, mean = 7.5, sd = 3.2)-dnorm(0, mean = 7.5, sd = 3.2) )/wl_de
            # },
            # constant = {
            #   wlc <- function(l) as.numeric((c1 <= l) & (l <= c2))^2
            #   wl_de <- sqrt(integrate(wlc, lower = 0, upper = maxLreal)$value)
            #   wl <- function(l) as.numeric((c1 <= l) & (l <= c2))/wl_de
            # }

    )
  } else {
    wlc <- function(l) wl.set(l)^2
    wl_de <- sqrt(integrate(wlc, lower = 0, upper = maxLreal)$value)
    wl <- function(l) wl.set(l)/wl_de
  }





  # model exposure process
  if(verbose) start = Sys.time()
  if(!interpolate) {
    SSx <- mgcv::smoothCon(s(t, bs = "bs", k = kx),
                           absorb.cons = FALSE,
                           data = data.frame(t = t))[[1]] ## reparameterize it later
    knots_x <- SSx$knots
    X <- SSx$X
    ## sum-to-zero reparameterization for SSx
    QRx <- qr(t(X) %*% as.vector(rep(1,nrow(X))))
    Qx <- qr.Q(QRx, complete = TRUE)
    Zx <- Qx[,2:ncol(Qx)]
    ## Check whether the design matrices are identical
    # X_repa <- SSx$X %*% Zx
    # max(unname(model.matrix(xt.fit))[,-1] - SSx$X %*% Zx) # SAME
  } else {
    # knots_x <- c(rep(0,2), t, rep(Nt+maxL+1,2))
    # X <- splines::splineDesign(knots = knots_x, x = c(0, t, Nt+maxL+1), outer.ok = TRUE)
    # X <- X[-c(1,Nt+2),]
    # knots_x <- c(rep(t[1]-1-0.2,3),t[1]-1-0.2, t[1]-1-0.001, t, Nt+maxL+0.5+0.001, Nt+maxL+0.5+0.2, rep(Nt+maxL+0.5+0.2,3))
    # knots_x <- c(rep(t[1]-1-0.2,4), t[1]-0.5-0.001, t, t[length(t)]+0.5+0.001, rep(t[length(t)]+1+0.2,4))
    knots_x <- c(t[1]-1-0.3, rep(t[1]-1-0.2,3), t[1]-0.5-0.001, t, t[length(t)]+0.5+0.001, rep(t[length(t)]+1+0.2,3), t[length(t)]+1+0.3)
    X <- splines::splineDesign(knots = knots_x, x = knots_x, outer.ok = TRUE)
    # X <- X[-c(1:2,Nt+maxL+7,Nt+maxL+8),]
  }




  if(!interpolate){
    xt.fit <- mgcv::gam(x~s(t, bs = "bs", k = kx), data = data.frame(x = PM25, t = t))
    ## coefficients for X(t)
    alpha_x <- xt.fit$coefficients

    # plot(t, xt.fit$fitted.values, type = "l", ylim = c(0,20))
    # lines(t, PM25, col = "red")
    Xt <- function(tnew) as.numeric(c(1, Bsplinevec2Con(tnew, knots_x, 4, Zx)) %*% alpha_x)
  } else {
    ## interpolate = TRUE
    ## set points for boundary and auxiliary boundary
    # alpha_x <- Interpolate(X, Zx, c(rep(0,4),PM25,rep(0,4)))
    Xsparse <- as(X, "dgCMatrix")
    # alpha_x <- Interpolate(Xsparse, c(rep(0,5),PM25,rep(0,5)))
    alpha_x1 <- Interpolate(Xsparse, c(rep(0,4), PM25[1], PM25, PM25[length(PM25)], rep(0,4)))
    xt.fit <- "interpolate"
    Xt1 <- function(tnew) as.numeric(Bsplinevec2(tnew, knots_x, 4) %*% alpha_x1)

    alpha_x2 <- Interpolate(Xsparse, c(rep(0,4), O3[1], O3, O3[length(O3)], rep(0,4)))
    Xt2 <- function(tnew) as.numeric(Bsplinevec2(tnew, knots_x, 4) %*% alpha_x2)

    alpha_x3 <- Interpolate(Xsparse, c(rep(0,4), NO2[1], NO2, NO2[length(NO2)], rep(0,4)))
    Xt3 <- function(tnew) as.numeric(Bsplinevec2(tnew, knots_x, 4) %*% alpha_x3)
  }
  if(verbose) cat("model exposure process: ", Sys.time() - start, "\n")




  ## The true exposure process is termed as the fitted function

  # CHECK: sapply(t, Xt) - xt.fit$fitted.values
  # Check: interpolate
  # Xtfit <- sapply(t, Xt)
  # max(abs(Xtfit - PM25)) #<1e-12



  ## True weighted exposure E = \int X(t-l) w(l) d l
  # E <- sapply(t, function(t.) integrate(Vectorize(function(l) wl(l)*Xt(t.-l)), lower = 0, upper = maxL)$value)
  l.eval <- seq(0, maxLreal, length.out = max(maxLreal, 500))
  wl.true <- sapply(l.eval, wl)

  wl.true.fit <- mgcv::gam(wl~s(l, bs = "bs", k = max(maxLreal, 100)), data = data.frame(wl = wl.true, l = l.eval))
  SSw.true <- mgcv::smoothCon(s(l, bs = "bs", k = max(maxLreal, 100)), absorb.cons = FALSE, data = data.frame(l = l.eval))[[1]]
  QRw.true <- qr(t(SSw.true$X) %*% as.vector(rep(1,nrow(SSw.true$X))))
  Qw.true <- qr.Q(QRw.true, complete = TRUE)
  Zw.true <- Qw.true[,2:ncol(Qw.true)]


  if(!interpolate) {
    if(verbose) start <- Sys.time()
    integ.E <- Integral(knots_x, SSw.true$knots, kx, max(maxLreal, 100), maxLreal, Zx, Zw.true, t.sim+0.5, alpha_x, TRUE)
    if(verbose) cat("integral 1: ", Sys.time() - start, "\n")
  }else {
    if(verbose) start <- Sys.time()
    integ.E1 <- acmeDLNM:::Integral_interpolate(knots_x, SSw.true$knots, kx, max(maxLreal, 100), maxLreal, Zw.true, t.sim+0.5, alpha_x1, TRUE)
    integ.E2 <- acmeDLNM:::Integral_interpolate(knots_x, SSw.true$knots, kx, max(maxLreal, 100), maxLreal, Zw.true, t.sim+0.5, alpha_x2, TRUE)
    integ.E3 <- acmeDLNM:::Integral_interpolate(knots_x, SSw.true$knots, kx, max(maxLreal, 100), maxLreal, Zw.true, t.sim+0.5, alpha_x3, TRUE)
    if(verbose) cat("integral 1: ", Sys.time() - start, "\n")
  }




  E.sim1 <- c(integ.E1$AlphaxD %*% wl.true.fit$coefficients)
  E.sim2 <- c(integ.E2$AlphaxD %*% wl.true.fit$coefficients)
  E.sim3 <- c(integ.E3$AlphaxD %*% wl.true.fit$coefficients)


  gamma <- gamma / sqrt(sum(gamma^2))
  # E.sim <- gamma[1] * E.sim1 + gamma[2] * E.sim2 + gamma[3] * E.sim3
  E.sim <- c(cbind(E.sim1, E.sim2, E.sim3) %*% gamma)

  ## association of weighted exposure f(E)
  Emax <- ceiling(max(E.sim))
  Emin <- floor(min(E.sim))


  if(is.null(fE.set)) {
    switch (fEtype,
            cubic = {
              fEtmp <- function(E) (E-(Emin + (Emax-Emin)*0.3))*(E-(Emin + (Emax-Emin)*0.2))*(E-(Emin + (Emax-Emin)*0.9))
              fE <- function(E) fEtmp(E) / (fEtmp(Emax)/2.5) + 2.5
            },
            linear = {
              fEtmp <- function(E) 0
              fE <- function(E) (E-(Emax+Emin)/2) / ((Emax - Emin)/3.5) + 3
            },
            # constant = {
            #   fEtmp <- function(E) 0
            #   fE <- function(E) 0
            # },
            quadratic = {
              fEtmp <- function(x){25*(dnorm(((2*((x-Emin)/(Emax-Emin)+0.18) - 1.1))))}
              fE <- function(x) -fEtmp(x)+11
            },
            simplelinear = {
              fEtmp <- function(E) 0
              fE <- function(E) 0.1*E
            }
            )
  } else {
    fEtmp <- function(E) 0
    fE <- function(E) fE.set(E)
  }

  ## linear predictor
  eta.sim <- sapply(E.sim, fE)
  Esurface.sim <- eta.sim - mean(eta.sim)
  if (!missingArg(other)) {
    eta.sim <- eta.sim + apply(other, 1, sum)
  }

  ## generate data
  y.sim <- sapply(eta.sim, function(eta.) rnbinom(n = 1, size = theta, mu = exp(eta.)))


  return(list(y.sim = y.sim,
              t.sim = t.sim,
              E.sim = E.sim,
              eta.sim = eta.sim,
              Esurface.sim = Esurface.sim,
              y = c(rep(0, maxL), y.sim),
              x1 = PM25,
              x1.orig = PM25.orig,
              x2 = O3,
              x2.orig = O3.orig,
              x3 = NO2,
              x3.orig = NO2.orig,
              gamma = gamma,
              t = t,
              xt.fit = xt.fit,
              alpha_x = list(alpha_x1, alpha_x2, alpha_x3), # for interpolation
              knots_x = knots_x, # for interpolation
              true.f = list(fE = fE, Emax = Emax, Emin = Emin, fEtmp = fEtmp,
                            wl = wl, wl_de = wl_de)
              ))
}
