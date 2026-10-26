#' Title Obtain predictions from the fitted model by acmeDLNM
#'
#' @param object object of class \code{acmeDLNM_fit}.
#' @param dat the data frame containing the variables.
#' @param CI.seed random seed for sampling, default \code{123}.
#' @param CI confidence level, default \code{0.95}.
#' @param ...
#' @details Predictions are on the log-mean scale for negative binomial fits
#'   and the mean scale for Gaussian fits. 
#'
#' @return
#' @importFrom mgcv s
#' @importFrom data.table as.data.table
#' @importFrom data.table setorder
#' @export
predict.acmeDLNM_fit <- function(object, dat, CI.seed = 123, CI = 0.95, ...) {

  family <- if(is.null(object$family)) "nb" else object$family
  gaussian <- identical(family, "gaussian")

  if(missingArg(dat)) {
    if (!is.null(object$eta)) {
      return(object$eta)
    } else {
      dat <- object$inputdata
    }
  }


  ## change character columns to factor. support random effects defined by mgcv::s(bs = "re")
  chr_col <- which(sapply(dat, class) == "character")
  if(length(chr_col) >= 1) {
    for(col. in chr_col) dat[,col.] <- factor(dat[,col.])
  }


  smooth = object$formula$smooth
  formula = object$formula$formula
  fe.cont = object$formula$fe.cont
  fe.varying = object$formula$fe.varying
  kw = object$data$kw
  kE = object$data$kE
  CI = object$CI
  maxL = object$data$maxL
  maxLreal = maxL+1





  sXobject <- eval(formula[[3]])

  sXdat <- dat

  ## change colnames in dataframe as x, t and y
  x.names.list <- sXobject$x

  colnames(sXdat)[which(colnames(sXdat) == sXobject$t)] <- "t"
  colnames(sXdat)[which(colnames(sXdat) == as.character(formula[[2]]))] <- "y"


  ## byvar
  byvar <- sXobject$by




  ### CONSTRUCTIONS #######
  ## time non-varying for group
  if(length(unique(sXdat$t)) < nrow(sXdat)) {

    if(is.null(fe.cont) & is.null(byvar)) stop("The exposure process data are duplicated at some time point. Please provide fe.cout.")

    if(!is.null(fe.cont)){
      # group_name.fe.cont <- fe.cont[[2]]
      # if(length(group_name.fe.cont) >= 2) {
      #   group_name.fe.cont <- as.character(group_name.fe.cont[2:length(group_name.fe.cont)])
      # } else {
      #   group_name.fe.cont <- as.character(group_name.fe.cont)
      # }
      group_name.fe.cont <- all.vars(fe.cont)
    } else {
      group_name.fe.cont <- NULL
    }
    if(!is.null(byvar)) {
      group_name.byvar <- byvar
    } else {
      group_name.byvar <- NULL
    }

    group_name <- c(group_name.fe.cont, group_name.byvar)

    ## split data
    sXdatlist <- split(as.data.table(sXdat), by = group_name) # use split in data.table to preserve order

  } else {
    sXdatlist <- list(sXdat)
  }

  ## set t starting at 1 and sort
  min_t <- min(sXdat$t)
  sXdatlist <- lapply(sXdatlist, function(sXdati) {
    sXdati$t <- sXdati$t - min_t + 1
    setorder(sXdati, t)
    return(data.frame(sXdati))
  })




  ## smooths from mgcv
  SSw <- object$smooth$wl

  SSfCon <- object$smooth$fE


  SwI_large <- object$data$SwI_large
  SwI <- object$data$SwI
  SfI <- object$data$SfI
  K <- object$data$K
  a <- object$data$a
  # Xfix <- object$data$Xfix
  knots_f <- object$data$knots_f

  knots_w <- object$data$knots_w
  Zf <- object$data$Zf
  Ufpen <- object$data$Ufpen
  Zfnew <- Zf %*% Ufpen
  Zw <- object$data$Zw
  Uwpen <- object$data$Uwpen
  Zwnew <- Zw %*% Uwpen
  kx.per500 <- object$data$kx.per500
  shift <- object$data$shift


  interpolate <- object$interpolate
  kx.per500 <- object$data$kx.per500


  ### 0. distributed lag term
  ## preparations for distributed lag terms
  DLprepare <- lapply(sXdatlist, function(sXdati){
    y <- sXdati$y
    t <- sXdati$t
    removed.t <- t[1:maxL] # removed id for the original t
    t <- t - min(t) + 1 # set t starting at 1

    DLi <- lapply(x.names.list, function(xx) {

      x <- sXdati[[xx]]
      
      Nti <- length(y) - maxL
      if((kx.per500 > 300) || (interpolate == TRUE)) {
        # if(verbose) cat("Interpolate the exposure process. \n")
        kx <- Nti+maxL + 4 + 2 # number of knots for X(t)
        # kx <- Nti+maxL + 4 # number of knots for X(t)
        interpolate <- TRUE
      } else{
        kx <- kx.per500 * ifelse(Nti < 500, 1, round(Nti/500)) # number of knots for X(t)
      }

      ### 0.1 Model exposure process
      # if(verbose) start <- Sys.time()
      if(!interpolate) {
        # if(verbose) cat("Start modelling exposure process ... ", "\n")
        SSx <- mgcv::smoothCon(s(t, bs = "bs", k = kx),
                               absorb.cons = FALSE,
                               data = data.frame(t = t))[[1]] ## reparameterize it later
        knots_x <- SSx$knots
        # X <- SSx$X

        ## sum-to-zero reparameterization for SSx
        QRx <- qr(t(SSx$X) %*% as.vector(rep(1,nrow(SSx$X))))
        Qx <- qr.Q(QRx, complete = TRUE)
        Zx <- Qx[,2:ncol(Qx)]
        ## Check whether the design matrices are identical
        # X_repa <- SSx$X %*% Zx
        # max(unname(model.matrix(xt.fit))[,-1] - SSx$X %*% Zx) # SAME

        xt.fit <- mgcv::gam(x~s(t, bs = "bs", k = kx), data = data.frame(x = x, t = t))
        ## coefficients for X(t)
        alpha_x <- xt.fit$coefficients
        # if(verbose) cat("modelling exposure process takes: ",
        #                 round(difftime(Sys.time(),start,units = 'secs'), 5),
        #                 "seconds. \n", "Consider setting interpolate = TRUE, if it (mgcv::gam function) takes too long. \n")
      } else {
        # knots_x <- c(rep(t[1]-1-0.2,4), t[1]-0.5-0.001, t, t[length(t)]+0.5+0.001, rep(t[length(t)]+1+0.2,4))
        knots_x <- c(t[1]-1-0.3, rep(t[1]-1-0.2,3), t[1]-0.5-0.001, t, t[length(t)]+0.5+0.001, rep(t[length(t)]+1+0.2,3), t[length(t)]+1+0.3)
        # knots_x <- c(rep(t[1]-1.5,3),t[1]-0.99, t, t[length(t)]+0.99, rep(t[length(t)]+1.5,3))
        X <- splines::splineDesign(knots = knots_x, x = knots_x, outer.ok = TRUE)
        ## interpolate = TRUE
        ## set points for boundary and auxiliary boundary
        Xsparse <- as(X, "dgCMatrix")
        alpha_x <- SIaceDLNM:::Interpolate(Xsparse, c(rep(0,4), x[1], x, x[length(x)], rep(0,4)))
        # alpha_x <- SIaceDLNM:::Interpolate(Xsparse, c(rep(0,4),x,rep(0,4)))
        xt.fit <- "interpolate"
        # if(verbose) cat("modelling exposure process takes: ",
        #                 round(difftime(Sys.time(),start,units = 'secs'), 5),
        #                 "seconds. \n")
      }


      ### 0.2 Integration
      t <- t[-(1:maxL)] # delete the first maxL days
      y <- y[-(1:maxL)] # delete the first maxL days
      x <- x[-(1:maxL)] # delete the first maxL days

      ### integration
      # if(verbose){
      #   cat("Start integration ... ", "\n")
      #   start <- Sys.time()
      # }
      if(!interpolate) {
        integral <- SIaceDLNM:::Integral(knots_x, knots_w, kx, kw, maxLreal, Zx, Zwnew, t+0.5, alpha_x, FALSE)
      } else {
        integral <- SIaceDLNM:::Integral_interpolate(knots_x, knots_w, kx, kw, maxLreal, Zwnew, t+0.5, alpha_x, FALSE)
      }
      # if(verbose) cat("integration takes: ",
      #                 round(difftime(Sys.time(),start,units = 'secs'), 5), "seconds. \n")

      ## linear predictor s(t) = f(\int w(l) X(t-l) dl) = f (B_inner(t) %*% alpha_w)
      ## where B_inner(t) = alpha_x %*% D, dim(D) = c(kx, kw), D_{p,q} = \int b_{xp}(t-l)b_{wq}(l) dl.
      ## b_{xp} is the p-th basis function for X(t)
      ## Stack B_inner(t) for t = (1:Nt)+40, we have B_inner. dim(B_inner) = c(Nt, kw)
      ## See Section X.X in Paper ...
      B_inner1 <- integral$AlphaxD

      ## For constraint 1: \int w(l)^2 dl = 1
      ## equivalent to w^T %*% Dw %*% w = 1
      Dw <- integral$Dw

      # E.maxi1 <- ceiling(max(integral$Xt2))
      E.maxi.square <- (integral$Xt2)^2

      return(list(B_inner = B_inner1,
                  Dw = Dw,
                  E.max.square = E.maxi.square))
    })

    return(list(B_inner = lapply(DLi, "[[", "B_inner"),
                Dw = DLi[[1]][["Dw"]],
                E.max = ceiling(max(sqrt(apply(do.call("cbind",lapply(DLi, "[[", "E.max.square")), 1, sum)))),
                removed.t =  removed.t
                ))
  })

  B_inner.list <- Reduce(function(x, y) Map(function(...){
      do.call(rbind, list(...))
    }, x, y), lapply(DLprepare, "[[", "B_inner"))



  Dw <- DLprepare[[1]]$Dw
  E.max <- object$data$E.max
  E.min <- object$data$E.min

  removed.t <- lapply(DLprepare, "[[", "removed.t")



  ## remove the starting time points
  sXdatlist <- mapply(function(sXdati, removed.ti) return(sXdati[-which(sXdati$t %in% removed.ti),]),
                      sXdatlist, removed.t, SIMPLIFY = FALSE)

  sXdat <- do.call("rbind", sXdatlist)


  # remove NA rows
  na.id <- is.na(sXdat$y)
  sXdat <- sXdat[!na.id, ]

  # B_inner <- list(B_inner1[!na.id, ], B_inner2[!na.id, ], B_inner3[!na.id, ])
  B_inner <- lapply(B_inner.list, function(xx) xx[!na.id, ])


  M <- length(B_inner) # number of exposures


  ## 1. time-varying fixed effects
  if(!is.null(fe.varying)) {
    Xlin <- stats::model.matrix(fe.varying,data=sXdat)[,-1,drop=F]
  } else {
    Xlin <- matrix(1, nrow = nrow(sXdat))
  }
  ## unpen.smooth = NULL
  # if(!is.null(unpen.smooth)) {
  #   unpen.smooth <- lapply(lapply(attr(terms(unpen.smooth),"term.labels"), function(text) parse(text = text)), eval)
  #   unpen.SS <- lapply(lapply(unpen.smooth,mgcv::smoothCon,data=sXdat,absorb.cons = TRUE),'[[',1) ## construct unpenalized smooth terms
  #   unpen.Xlist <- lapply(unpen.SS,'[[','X')
  #   unpen.X <- Reduce(cbind,unpen.Xlist) ## design matrix
  #   Xlin <- cbind(Xlin, unpen.X)
  # }
  ### code following https://github.com/awstringer1/mam/blob/master/R/mam.R
  ## start following ##
  ## 2.2 smooth term
  numsmooth <- 0 # initialize
  if(!is.null(smooth)){
    smooth <- lapply(lapply(attr(terms(smooth),"term.labels"), function(text) parse(text = text)), eval)

    SS <- lapply(lapply(smooth,mgcv::smoothCon,data=sXdat,absorb.cons = TRUE),'[[',1) ## construct smooth terms
    numsmooth <- length(smooth) # Number of smooth terms

    EE <- lapply(lapply(lapply(SS,'[[','S'),'[[',1),eigen) ## eigen decomposition for penalty matrix

    p <- sapply(lapply(EE,'[[','vectors'),ncol) ## dimension of penalty matrix
    r <- sapply(lapply(EE,'[[','values'),function(x) sum(x>.Machine$double.eps)) ## rank of penalty matrix
    m <- p-r ## dim of null space (minus intercept)
    URlist <- mapply(function(x,y) x[ ,1:y],lapply(EE,'[[','vectors'),r,SIMPLIFY = FALSE)
    UFlist <- mapply(
      function(x,y,z) {
        if (y<z) return(x[ ,(1+y):z])
        newsparsemat(z,z)
      },lapply(EE,'[[','vectors'),r,p,SIMPLIFY = FALSE)
    URlist <- lapply(URlist,cbind) # Ensure they stay matrices
    UFlist <- lapply(UFlist,cbind) # Ensure they stay matrices

    UR <- Matrix::bdiag(URlist)
    UF <- Matrix::bdiag(UFlist)
    # if m=1 UF gets coerced to numeric
    if (!is.matrix(UF)) UF <- cbind(UF)

    Dpi <- Matrix::Diagonal(sum(r),1 / sqrt(Reduce(c,lapply(lapply(EE,'[[','values'),function(x) x[x>.Machine$double.eps]))))

    Xlist <- lapply(SS,'[[','X')
    X <- Reduce(cbind,Xlist) ## design matrix

    Xrand <- as.matrix(X %*% UR %*% Dpi) ## reparametrized
    Xfix <- as.matrix(X %*% UF)
    dups <- !base::duplicated(t(Xfix)) & apply(Xfix,2,function(x) !all(x==0)) # Remove the duplicated intercepts
    if (length(dups) > 1) Xfix <- Xfix[ ,which(dups)]

    model.choice = "with.smooth"
  } else{
    model.choice = "without.smooth"
    Xrand <- matrix(numeric(0), nrow = nrow(sXdat), ncol = 0)
    Xfix <- matrix(numeric(0), nrow = nrow(sXdat), ncol = 0)
    r <- numeric(0)
  }


  # 2.3 add the intercept
  if(exists("Xfix")) {
    Xfix <- cbind(Xlin,Xfix) ## linear effect + unpenalized columns
  } else {
    Xfix <- Xlin
  }
  if (any(apply(Xfix,2,function(x) all(x==1)))) Xfix <- Xfix[,-(apply(Xfix,2,function(x) all(x==1)))]

  if (!is.null(fe.cont)){
    Xgroup <- stats::model.matrix(fe.cont,data=sXdat)[,-1,drop=FALSE]
    if (any(apply(Xgroup,2,function(x) all(x==1)))) Xgroup <- Xgroup[,-(apply(Xgroup,2,function(x) all(x==1)))] # remove the intercept
    Xfix <- cbind(1, Xgroup, Xfix) # add the intercept to the first column
  } else {
    Xfix <- cbind(1,Xfix) # add the intercept to the first column
  }


  ### END following https://github.com/awstringer1/mam/blob/master/R/mam.R #######

  ## 2.4 offset
  ## TODO: update offset
  # if(!is.null(offset)) {
  #   Xoffset.original <- as.vector(stats::model.matrix(offset,data=sXdat)[,-1,drop=F])
  #   Xoffset.original.min <- min(Xoffset.original)
  #   Xoffset <- Xoffset.original - Xoffset.original.min
  # } else {
  #   Xoffset.original <- as.vector(rep(0, nrow(sXdat)))
  #   Xoffset.original.min <- 0
  #   Xoffset <- as.vector(rep(0, nrow(sXdat)))
  # }
  Xoffset <- if(gaussian) rep(0, nrow(sXdat)) else object$data$offset$Xoffset


  y <- sXdat$y


  B_inner_sum <- Reduce("+", mapply("*", B_inner, object$point$index, SIMPLIFY = FALSE))
  E <- B_inner_sum %*% object$point$alpha_w

  Bf <- sapply(E, function(Ei) Bsplinevec2Con(Ei, knots_f, 4, Zfnew))
  eta1 <- as.vector(t(Bf) %*% object$point$alpha_f)
  eta2 <- as.vector(Xfix %*% object$point$betaF)
  if(exists("Xrand")) eta2 <- eta2 + as.vector(Xrand %*% object$point$betaR)

  eta.est <- eta1 + eta2 + Xoffset


  eta_E.est <- eta1 + mean(eta2)
  eta_other.est <- eta2 - mean(eta2)



  Bindex <- object$data$Bindex
  if(is.null(Bindex)) {
    cindex <- rep(1, M)
    Iindex <- diag(1, nrow = M, ncol = M) # the so-called basis function for index weight
    Vindex <- Iindex %*% cindex
    Qfindex <- qr.Q(qr(Vindex), complete = TRUE)
    Zfindex <- Qfindex[,2:ncol(Qfindex)]
    Bindex <- cbind(cindex, Iindex %*% Zfindex)
  }


  # load("test.RData")
  # acmeDLNMbuild <- acmeDLNM:::acmeDLNMbuild
  # acmeDLNMCI <- acmeDLNM:::acmeDLNMCI

  mod.build <- acmeDLNMbuild(y, B_inner, knots_f, SwI_large, SfI, Dw,
                            Xrand, Xfix, Zfnew, Xoffset, r,
                            K,a,
                            Bindex,
                            object$point$alpha_f,
                            object$point$phi,
                            object$point$delta,
                            object$point$log_theta, object$point$log_smoothing_f, object$point$log_smoothing_w,
                            object$point$betaR, object$point$betaF, object$point$log_smoothing,
                            gaussian = gaussian)

  mod.address <- mod.build$address.eigen

  R.CI <- nrow(object$CI.sample[[1]])

  # R.CI <- 2000
  sampled <- acmeDLNMCI(mod.address,
                        R.CI, CI.seed,
                        TRUE, FALSE, FALSE,
                        object$Hessian_inner)
  eta.CI <- apply(sampled$eta_sample_mat, 2, function(row.) quantile(row., c((1-CI)/2, 1-(1-CI)/2)))
  out <- list(eta = data.frame(est = eta.est,
                        ll = eta.CI[1,],
                        ul = eta.CI[2,]),
              sample = sampled$eta_sample_mat)

  return(out)
}
