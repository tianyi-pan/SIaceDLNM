#' Title
#'
#' @param delta
#' @param formula.mgcv
#' @param sXdat
#' @param opt
#' @param verbose
#' @param family outcome family, \code{"nb"} (default) or \code{"gaussian"}.
#'
#' @returns
#' @importFrom mgcv bam
#' @importFrom mgcv s
#' @importFrom mgcv nb
#'
#' @examples
si <- function(delta, formula.mgcv, a.names.vec, sXdat, Bindex, opt = TRUE, verbose, family = "nb") {
  if(!is.character(family) || length(family) != 1L || is.na(family) ||
     !family %in% c("nb", "gaussian")) stop("family must be 'nb' or 'gaussian'.")
  alpha <- c(1,delta)
  # kk <- sqrt(sum(alpha^2))
  kk <- as.numeric(t(alpha) %*% t(Bindex) %*% Bindex %*% alpha)
  alpha <- Bindex %*% alpha/kk  ## so now ||alpha||=1
  a <- as.vector(as.matrix(sXdat[a.names.vec]) %*% alpha)
  # a <- alpha * sXdat$x1.avglag + alpha[2] * sXdat$x2.avglag + alpha[3] * sXdat$x3.avglag

  sXdat$a <- a

  b.opt <- bam(formula.mgcv,
               data = sXdat,
               family=if(identical(family, "gaussian")) stats::gaussian(link = "identity") else nb(),
               discrete = TRUE) ## fit model
  if(verbose) cat(".")
  if (opt)
    return(b.opt$gcv.ubre)
  else {
    b.opt$alpha <- alpha  ## add alpha
    # J <- outer(alpha,-delta/kk^2) ## compute Jacobian
    # for (j in 1:length(delta)) J[j+1,j] <- J[j+1,j] + 1/kk
    # b.opt$J <- J ## dalpha_i/ddelta_j
    return(b.opt)
  }
}


#' Title
#'
#' @param formula.list
#' @param kE
#' @param sXdat
#' @param M
#' @param verbose
#' @param family outcome family, \code{"nb"} (default) or \code{"gaussian"}.
#'
#' @returns
#' @export
#'
#' @examples
sim <- function(formula.list, x.names.list, kE, sXdat, Bindex, M, maxL, verbose, con.index.init, family = "nb") {

  if(!is.character(family) || length(family) != 1L || is.na(family) ||
     !family %in% c("nb", "gaussian")) stop("family must be 'nb' or 'gaussian'.")

  wl.discrete <- 1 - (0:maxL)/maxL

  wl.discrete <- wl.discrete/sqrt(sum(wl.discrete^2))


  a.names.vec <- paste0("a", seq_len(length(x.names.list)))

  for(i in seq_len(length(x.names.list))){
    sXdat[[a.names.vec[i]]] <- apply(sapply(1:(maxL+1), function(ii) {
      dplyr::lag(sXdat[[x.names.list[[i]]]], ii-1) * wl.discrete[ii]
    }), 1, sum)
  }
#
#
#   sXdat$x2.avglag <- apply(sapply(1:(maxL+1), function(ii) {
#     dplyr::lag(sXdat$x2, ii-1) * wl.discrete[ii]
#   }), 1, sum)
#
#   sXdat$x3.avglag <- apply(sapply(1:(maxL+1), function(ii) {
#     dplyr::lag(sXdat$x3, ii-1) * wl.discrete[ii]
#   }), 1, sum)



  if(verbose) cat("Fitting a single index model using mgcv::bam() to obtain init ... \n")
  formula.list.mgcv <- Filter(Negate(is.null), formula.list[names(formula.list) %in% c("smooth","fe.varying")])

  formula.mgcv.unpen <- paste0("y~s(a, bs='bs', fx = TRUE, k = ", kE, ")")
  formula.mgcv.pen <- paste0("y~s(a, bs='bs', fx = FALSE, k = ", kE, ")")

  if(length(formula.list.mgcv) > 0) {
    formula.other.mgcv <- paste(formula.list.mgcv, collapse = "+")
    formula.other.mgcv <- gsub("~", "", formula.other.mgcv)
    formula.other.mgcv <- gsub("\"", "'", formula.other.mgcv)
    formula.mgcv.unpen <- paste0(formula.mgcv.unpen, "+", formula.other.mgcv)
    formula.mgcv.pen <- paste0(formula.mgcv.pen, "+", formula.other.mgcv)
  }
  formula.mgcv.unpen <- as.formula(formula.mgcv.unpen)
  formula.mgcv.pen <- as.formula(formula.mgcv.pen)


  if(M == 2) {
    method = "Brent"
    lower = -10
    upper = 10
  } else {
    method = "Nelder-Mead"
    lower = -Inf
    upper = Inf
  }
  if(missingArg(con.index.init)) {
    th0 <- rep(1,M-1)
    f0 <- optim(th0, si, formula.mgcv = formula.mgcv.unpen,
                a.names.vec = a.names.vec,
                sXdat = sXdat, Bindex = Bindex, verbose = verbose, family = family,
                method = method, lower = lower, upper = upper)



    f1 <- optim(f0$par, si, formula.mgcv = formula.mgcv.pen,
                a.names.vec = a.names.vec,
                sXdat = sXdat, Bindex = Bindex, verbose = verbose, family = family,
                method = method, lower = lower, upper = upper)


    delta.est <- f1$par

    b.opt <- si(delta = delta.est, formula.mgcv = formula.mgcv.pen,
                a.names.vec = a.names.vec,
                sXdat = sXdat, Bindex = Bindex, opt = FALSE, verbose = verbose, family = family)
  } else {
    delta.est <- con.index.init
    b.opt <- si(delta = delta.est, formula.mgcv = formula.mgcv.pen,
                a.names.vec = a.names.vec,
                sXdat = sXdat, Bindex = Bindex, opt = FALSE, verbose = verbose, family = family)
  }

  if(verbose) cat("Single index model: Done! \n")

  return(list(delta.init = delta.est,
              mod.mgcv = b.opt)
         )
}


### TODO: add Bindex in the code


