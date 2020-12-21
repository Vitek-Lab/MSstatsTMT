#############################################################################
#    This file is revised from the lmerTest package for R
#   https://github.com/rforge/lmertest/blob/master/pkg/lmerTest/R
#############################################################################

##########################################################################
## Create rho vector containing info about mixed model ##################
##########################################################################
#' @importFrom lme4 fixef getME
#' @importFrom stats sigma
#' @keywords internal
.rhoInit <- function(rho, model, change.contr = FALSE, mf.final = NULL)
{
  if(change.contr)
    rho$model <- .updateModel(model, mf.final = mf.final, change.contr = change.contr)
  else
    rho$model <- model
  rho$fixEffs <- lme4::fixef(rho$model)
  rho$sigma <- stats::sigma(rho$model)
  rho$thopt <- lme4::getME(rho$model, "theta")
  return(rho)  
}

#############################################
## returns Lc %*% vcov as a function of theta parameters %*% t(Lc)
#############################################
#' @importFrom lme4 isGLMM isLMM fixef
#' @importFrom methods is
#' @keywords internal
.vcovLThetaL <- function(fm) {
  stopifnot(is(fm, "merMod"))
  
  np <- length(fm@pp$theta)
  nf <- length(lme4::fixef(fm))
  if (!lme4::isGLMM(fm)) 
    np <- np + 1L

  ff2 <- .updateModel(fm, devFunOnly = TRUE) 
  
  envff2 <- environment(ff2)
  
  if (lme4::isLMM(fm)) {
    ans <- function(Lc, thpars) {
      stopifnot(is.numeric(thpars), length(thpars) == np)
      
      sigma2 <- thpars[np]^2
      ff2(thpars[-np])
      
      vcov_unscaled <- tcrossprod(envff2$pp$RXi()) 
      vcov_out <- sigma2 * vcov_unscaled
      
      return(list(varcor = as.matrix(Lc %*% as.matrix(vcov_out) %*% t(Lc)),
                  unscaled.varcor = vcov_unscaled,
                  sigma2 = sigma2)) 
    }
  } 
  class(ans) <- ".vcovLThetaL"
  
  ans
}

#############################################
## calculate gradient
#############################################
#' @keywords internal
.mygrad <- function(fun, x, delta = 1e-4,
                   method = c("central", "forward", "backward"), ...)
{
  method <- match.arg(method)
  nx <- length(x)
  if(method %in% c("central", "forward")) {
    Xadd <- matrix(rep(x, nx), nrow=nx, byrow=TRUE) + diag(delta, nx)
    fadd <- apply(Xadd, 1, fun, ...)
  }
  if(method %in% c("central", "backward")) {
    Xsub <- matrix(rep(x, nx), nrow=nx, byrow=TRUE) - diag(delta, nx)
    fsub <- apply(Xsub, 1, fun, ...) ## eval.parent perhaps?
  }
  res <- switch(method,
                "forward" = (fadd - fun(x, ...)) / delta,
                "backward" = (fun(x, ...) - fsub) / delta,
                "central" = (fadd - fsub) / (2 * delta)
  )
  res
}

#############################################
## devfun function as a function of optimal parameters
#############################################
#' @importFrom stats formula getCall terms update.formula
#' @keywords internal
.updateModel <- function(object, mf.final = NULL, ..., change.contr = FALSE) {
  if (is.null(call <- getCall(object)))
    stop("object should contain a 'call' component")
  extras <- match.call(expand.dots = FALSE)$...
  if (!is.null(mf.final))
    call$formula <- update.formula(formula(object), mf.final)
  if(any(grepl("sample", call))){
    call <- as.list(call)[-which(names(as.list(call)) %in% c("data", "subset"))]
    call[["data"]] <- quote(model.frame(object))
  }
  if (length(extras) > 0) {
    existing <- !is.na(match(names(extras), names(call)))
    for (a in names(extras)[existing]) call[[a]] <- extras[[a]]
    if (any(!existing)) {
      call <- c(as.list(call), extras[!existing])
    }
  }
  if(change.contr){
    mm <- model.matrix(object)
    contr <- attr(mm,"contrasts")
    ## change contrasts for F tests calculations
    ## list of contrasts for factors
    if(change.contr && length(which(unlist(contr)!="contr.SAS")) > 0)
    {
      names.facs <- names(contr)
      l.lmerTest.private.contrast <- as.list(rep("contr.SAS",length(names.facs)))
      names(l.lmerTest.private.contrast) <- names(contr)
      call[["contrasts"]] <- l.lmerTest.private.contrast
    }
    else if(!is.null(contr)){
      call[["contrasts"]] <- contr[names(contr) %in% attr(terms(call$formula),"term.labels")]
    }
  }
  call <- as.call(call)
  ff <- environment(formula(object))
  pf <- parent.frame()  ## save parent frame in case we need it
  sf <- sys.frames()[[1]]
  ff2 <- environment(object)
  tryCatch(eval(call, envir=ff),
           error=function(e) {
             tryCatch(eval(call, envir=sf),
                      error=function(e) {
                        tryCatch(eval(call, envir=pf),
                                 error=function(e) {
                                   eval(call, envir=ff2)
                                 })})})
}

#' # retired fuction (2020.12.17)#############################################
#' ## calculates asymptotic variance covariance matrix of variance parameters based on theta
#' #############################################
#' #' @keywords internal
#' .calcApvar <- function(rho){
#'   ## based on theta parameters and sigma
#'   dd <- .devfunTheta(rho$model)
#'   h <- .myhess(dd, c(rho$thopt, sigma = rho$sigma))
#' 
#'   ch <- try(chol(h), silent=TRUE)
#'   if(inherits(ch, "try-error")) {
#'     return(rho)
#'   }
#'   A <- 2*chol2inv(ch)
#' 
#'   eigval <- eigen(h, symmetric=TRUE, only.values=TRUE)$values
#'   isposA <- TRUE
#'   if(min(eigval) < sqrt(.Machine$double.eps)) ## tol ~ sqrt(.Machine$double.eps)
#'     isposA <- FALSE
#' 
#'   if(!isposA)
#'     print("Asymptotic covariance matrix A is not positive!")
#' 
#'   A
#' }
#' 
#' #############################################
#' ## devfun function as a function of optimal parameters
#' #############################################
#' #' @importFrom lme4 isGLMM isLMM getME
#' #' @importFrom methods is
#' #' @keywords internal
#' .devfunTheta <- function(fm)
#' {
#'   stopifnot(is(fm, "merMod"))
#' 
#'   np <- length(fm@pp$theta)
#'   nf <- length(fixef(fm))
#'   if (!lme4::isGLMM(fm))
#'     np <- np + 1L
#'   n <- nrow(fm@pp$V)
#' 
#'   ff <- .updateModel(fm, devFunOnly = TRUE)
#'   reml <- lme4::getME(fm, "is_REML")
#' 
#'   envff <- environment(ff)
#' 
#'   if (lme4::isLMM(fm)) {
#'     ans <- function(thpars) {
#'       stopifnot(is.numeric(thpars), length(thpars) == np)
#' 
#'       ff(thpars[-np])
#' 
#'       sigsq <- thpars[np]^2
#'       dev <- envff$pp$ldL2() + (envff$resp$wrss() + envff$pp$sqrL(1))/sigsq + n * log(2 * pi * sigsq)
#'       if(reml){
#'         p <- ncol(envff$pp$RX())
#'         dev <- dev + 2*determinant(envff$pp$RX())$modulus - p * log(2 * pi * sigsq)
#'       }
#'       return(dev)
#'     }
#'   }
#' 
#'   attr(ans, "thopt") <- fm@pp$theta
#'   class(ans) <- ".devfunTheta"
#'   ans
#' }
#' 
#' #############################################
#' ## calculate hessian matrix
#' #############################################
#' #' @keywords internal
#' .myhess <- function(fun, x, fx=NULL, delta=1e-4, ...) {
#'   nx <- length(x)
#'   fx <- if(!is.null(fx)) fx else fun(x, ...)
#'   H <- array(NA, dim=c(nx, nx))
#'   for(j in 1:nx) {
#'     ## Diagonal elements:
#'     xadd <- xsub <- x
#'     xadd[j] <- x[j] + delta
#'     xsub[j] <- x[j] - delta
#'     H[j, j] <- (fun(xadd, ...) - 2 * fx +
#'                   fun(xsub, ...)) / delta^2
#'     ## Upper triangular (off diagonal) elements:
#'     for(i in 1:nx) {
#'       if(i >= j) break
#'       xaa <- xas <- xsa <- xss <- x
#'       xaa[c(i, j)] <- x[c(i, j)] + c(delta, delta)
#'       xas[c(i, j)] <- x[c(i, j)] + c(delta, -delta)
#'       xsa[c(i, j)] <- x[c(i, j)] + c(-delta, delta)
#'       xss[c(i, j)] <- x[c(i, j)] - c(delta, delta)
#'       H[i, j] <- (fun(xaa, ...) - fun(xas, ...) -
#'                     fun(xsa, ...) + fun(xss, ...)) /
#'         (4 * delta^2)
#'     }
#'   }
#'   ## Fill in lower triangle:
#'   H[lower.tri(H)] <- t(H)[lower.tri(H)]
#' 
#'   H
#' }