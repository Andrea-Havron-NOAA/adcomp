## Copyright (C) 2013-2015 Kasper Kristensen
## License: GPL-2

##' After optimization of an AD model, \code{sdreport} is used to
##' calculate standard deviations of all model parameters, including
##' non linear functions of random effects and parameters specified
##' through the ADREPORT() macro from the user template.
##'
##' First, the Hessian wrt. the parameter vector (\eqn{\theta}) is
##' calculated.  The parameter covariance matrix is approximated by
##' \deqn{V(\hat\theta)=-\nabla^2 l(\hat\theta)^{-1}} where \eqn{l}
##' denotes the log likelihood function (i.e. \code{-obj$fn}).  If
##' \code{ignore.parm.uncertainty=TRUE} then the Hessian calculation
##' is omitted and a zero-matrix is used in place of
##' \eqn{V(\hat\theta)}.
##'
##' For non-random effect models the standard delta-method is used to
##' calculate the covariance matrix of transformed parameters. Let
##' \eqn{\phi(\theta)} denote some non-linear function of
##' \eqn{\theta}. Then \deqn{V(\phi(\hat\theta))\approx \nabla\phi
##' V(\hat\theta) \nabla\phi'}
##'
##' The covariance matrix of reported variables
##' \eqn{V(\phi(\hat\theta))} is returned by default. This can cause
##' high memory usage if many variables are ADREPORTed. Use
##' \code{getReportCovariance=FALSE} to only return standard errors.
##' In case standard deviations are not required one can completely skip
##' the delta method using \code{skip.delta.method=TRUE}.
##'
##' For random effect models a generalized delta-method is used. First
##' the joint covariance of random effect and parameter estimation error is approximated
##' by
##' \deqn{V \left( \begin{array}{cc} \hat u - u \cr \hat\theta - \theta \end{array} \right) \approx
##' \left( \begin{array}{cc} H_{uu}^{-1} & 0 \cr 0 & 0 \end{array} \right) +
##' J V(\hat\theta) J'
##' }
##' where \eqn{H_{uu}} denotes random effect block of the full joint
##' Hessian of \code{obj$env$f} and \eqn{J} denotes the Jacobian of
##' \eqn{\left( \begin{array}{cc}\hat u(\theta) \cr \theta \end{array} \right)} wrt. \eqn{\theta}.
##' Here, the first term represents the expected conditional variance
##' of the estimation error given the data and the second term represents the variance
##' of the conditional mean of the estimation error given the data.
##'
##' Now the delta method can be applied on a general non-linear
##' function \eqn{\phi(u,\theta)} of random effects \eqn{u} and
##' parameters \eqn{\theta}:
##' \deqn{V\left(\phi(\hat u,\hat\theta) - \phi(u,\theta) \right)\approx \nabla\phi V \left( \begin{array}{cc}
##' \hat u - u \cr \hat\theta - \theta \end{array} \right) \nabla\phi'}
##'
##' The full joint covariance is not returned by default, because it
##' may require large amounts of memory.  It may be obtained by
##' specifying \code{getJointPrecision=TRUE}, in which case \eqn{V
##' \left( \begin{array}{cc} \hat u - u \cr \hat\theta - \theta \end{array} \right) ^{-1} } will be part of the
##' output. This matrix must be manually inverted using
##' \code{solve(jointPrecision)} in order to get the joint covariance
##' matrix. Note, that the parameter order will follow the original
##' order (i.e. \code{obj$env$par}).
##'
##' Using \eqn{\phi(\hat u,\theta)} as estimator of
##' \eqn{\phi(u,\theta)} may result in substantial bias. This may be
##' the case if either \eqn{\phi} is non-linear or if the distribution
##' of \eqn{u} given \eqn{x} (data) is sufficiently non-symmetric.  A
##' generic correction is enabled with \code{bias.correct=TRUE}. It is
##' based on the identity
##' \deqn{E_{\theta}[\phi(u,\theta)|x] =
##' \partial_\varepsilon\left(\log \int \exp(-f(u,\theta) +
##' \varepsilon \phi(u,\theta))\:du\right)_{|\varepsilon=0}}
##' stating that the conditional expectation can be written as a
##' marginal likelihood gradient wrt. a nuisance parameter
##' \eqn{\varepsilon}.
##' The marginal likelihood is replaced by its Laplace approximation.
##'
##' If \code{bias.correct.control$sd=TRUE} the variance of the
##' estimator is calculated using
##' \deqn{V_{\theta}[\phi(u,\theta)|x] =
##' \partial_\varepsilon^2\left(\log \int \exp(-f(u,\theta) +
##' \varepsilon \phi(u,\theta))\:du\right)_{|\varepsilon=0}}
##' A further correction is added to this variance to account for the
##' effect of replacing \eqn{\theta} by the MLE \eqn{\hat\theta}
##' (unless \code{ignore.parm.uncertainty=TRUE}).
##'
##' Bias correction can be be performed in chunks in order to reduce
##' memory usage or in order to only bias correct a subset of
##' variables. First option is to pass a list of indices as
##' \code{bias.correct.control$split}. E.g. a list
##' \code{list(1:2,3:4)} calculates the first four ADREPORTed
##' variables in two chunks.
##' The internal function \code{obj$env$ADreportIndex()}
##' gives an overview of the possible indices of ADREPORTed variables.
##'
##' Second option is to pass the number of
##' chunks as \code{bias.correct.control$nsplit} in which case all
##' ADREPORTed variables are bias corrected in the specified number of
##' chunks.
##' Also note that \code{skip.delta.method} may be necessary when bias
##' correcting a large number of variables.
##'
##' @title General sdreport function.
##' @param obj Object returned by \code{MakeADFun}
##' @param par.fixed Optional. Parameter estimate (will be known to \code{obj} when an optimization has been carried out).
##' @param hessian.fixed Optional. Hessian wrt. parameters (will be calculated from \code{obj} if missing).
##' @param getJointPrecision Optional. Return full joint precision matrix of random effects and parameters?
##' @param bias.correct logical indicating if bias correction should be applied
##' @param bias.correct.control a \code{list} of bias correction options; currently \code{sd}, \code{split} and \code{nsplit} are used - see details.
##' @param ignore.parm.uncertainty Optional. Ignore estimation variance of parameters?
##' @param getReportCovariance Get full covariance matrix of ADREPORTed variables?
##' @param skip.delta.method Skip the delta method? (\code{FALSE} by default)
##' @return Object of class \code{sdreport}
##' @seealso \code{\link{summary.sdreport}}, \code{\link{print.sdreport}}, \code{\link{as.list.sdreport}}
##' @examples
##' \dontrun{
##' runExample("linreg_parallel", thisR = TRUE) ## Non-random effect example
##' sdreport(obj) }
##'
##' runExample("simple", thisR = TRUE)          ## Random effect example
##' rep <- sdreport(obj)
##' summary(rep, "random")                      ## Only random effects
##' summary(rep, "fixed", p.value = TRUE)       ## Only non-random effects
##' summary(rep, "report")                      ## Only report
##'
##' ## Bias correction
##' rep <- sdreport(obj, bias.correct = TRUE)
##' summary(rep, "report")                      ## Include bias correction
sdreport <- function(obj,par.fixed=NULL,hessian.fixed=NULL,getJointPrecision=FALSE,bias.correct=FALSE,
                     bias.correct.control=list(sd=FALSE, split=NULL, nsplit=NULL), ignore.parm.uncertainty = FALSE,
                     getReportCovariance=TRUE, skip.delta.method=FALSE){
  if(is.null(obj$env$ADGrad) & (!is.null(obj$env$random)))
    stop("Cannot calculate sd's without type ADGrad available in object for random effect models.")
  ## Make object to calculate ADREPORT vector
  obj2 <- MakeADFun(obj$env$data,
                    obj$env$parameters,
                    type = "ADFun",
                    ADreport = TRUE,
                    DLL = obj$env$DLL,
                    silent = obj$env$silent)
  r <- obj$env$random
  ## Get full parameter (par), Fixed effects parameter (par.fixed)
  ## and fixed effect gradient (gradient.fixed)
  if(is.null(par.fixed)){ ## Parameter estimate not specified - use best encountered parameter
    par <- obj$env$last.par.best
    if(!is.null(r))par.fixed <- par[-r] else par.fixed <- par
    gradient.fixed <- obj$gr(par.fixed)
  } else {
    gradient.fixed <- obj$gr(par.fixed) ## <-- updates last.par
    par <- obj$env$last.par
  }
  ## In case of empty parameter vector:
  if(length(par.fixed)==0) ignore.parm.uncertainty <- TRUE
  ## Get Hessian wrt. fixed effects (hessian.fixed) and check if positive definite (pdHess).
  if(ignore.parm.uncertainty){
      hessian.fixed <- NULL
      pdHess <- TRUE
      Vtheta <- matrix(0, length(par.fixed), length(par.fixed))
  } else {
      if(is.null(hessian.fixed)){
          hessian.fixed <- optimHess(par.fixed,obj$fn,obj$gr) ## Marginal precision of theta.
      }
      pdHess <- !is.character(try(chol(hessian.fixed),silent=TRUE))
      Vtheta <- try(solve(hessian.fixed),silent=TRUE)
      if(is(Vtheta, "try-error")) Vtheta <- hessian.fixed * NaN
  }
  ## Get random effect block of the full joint Hessian (hessian.random) and its
  ## Cholesky factor (L)
  if(!is.null(r)){
    hessian.random <- obj$env$spHess(par,random=TRUE)   ## Conditional prec. of u|theta
    L <- obj$env$L.created.by.newton
    if(!is.null(L)){ ## Re-use symbolic factorization if exists
      updateCholesky(L,hessian.random)
      hessian.random@factors <- list(SPdCholesky=L)
    }
  }
  ## Get ADreport vector (phi)
  phi <- try(obj2$fn(par), silent=TRUE)    ## NOTE_1: obj2 forward sweep now initialized !
  if(is.character(phi) | length(phi)==0){
      phi <- numeric(0)
  }
  ADGradForward0Initialized <- FALSE
  ADGradForward0Initialize <- function() { ## NOTE_2: ADGrad forward sweep now initialized !
      obj$env$f(par, order = 0, type = "ADGrad")
      ADGradForward0Initialized <<- TRUE
  }
  doDeltaMethod <- function(chunk=NULL){
      ## ======== Determine case
      ## If no random effects use standard delta method
      simpleCase <- is.null(r)
      if(length(phi)==0){ ## Nothing to report
          simpleCase <- TRUE
      } else { ## Something to report - get derivatives
          if(is.null(chunk)){ ## Do all at once
              Dphi <- obj2$gr(par)
          } else {
              ## Do *chunk* only
              ## Reduce to Dphi[chunk,] and phi[chunk]
              w <- rep(0, length(phi))
              phiDeriv <- function(i){
                  w[i] <- 1
                  obj2$env$f(par, order=1, rangeweight=w, doforward=0) ## See NOTE_1
              }
              Dphi <- t( sapply(chunk, phiDeriv) )
              phi <- phi[chunk]
          }
          if(!is.null(r)){
              Dphi.random <- Dphi[,r,drop=FALSE]
              Dphi.fixed <- Dphi[,-r,drop=FALSE]
              if(all(Dphi.random==0)){ ## Fall back to simple case
                  simpleCase <- TRUE
                  Dphi <- Dphi.fixed
              }
          }
      }
      ## ======== Do delta method
      ## Get covariance (cov)
      if(simpleCase){
          if(length(phi)>0){
              cov <- Dphi %*% Vtheta %*% t(Dphi)
          } else cov <- matrix(,0,0)
      } else {
          tmp <- solve(hessian.random,t(Dphi.random))
          tmp <- as.matrix(tmp)
          term1 <- Dphi.random%*%tmp ## first term.
          if(ignore.parm.uncertainty){
              term2 <- 0
          } else {
              ## Use columns of tmp as direction for reverse mode sweep
              f <- obj$env$f
              w <- rep(0, length(par))
              if(!ADGradForward0Initialized) ADGradForward0Initialize()
              reverse.sweep <- function(i){
                  w[r] <- tmp[,i]
                  -f(par, order = 1, type = "ADGrad", rangeweight = w, doforward=0)[-r]
              }
              A <- t(do.call("cbind",lapply(seq_along(phi), reverse.sweep))) + Dphi.fixed
              term2 <- A %*% (Vtheta %*% t(A)) ## second term
          }
          cov <- term1 + term2
      }
      ##list(phi=phi, cov=cov)
      cov
  }
  if (!skip.delta.method) {
      if (getReportCovariance) { ## Get all
          cov <- doDeltaMethod()
          sd <- sqrt(diag(cov))
      } else {
          tmp <- lapply(seq_along(phi), doDeltaMethod)
          sd <- sqrt(as.numeric(unlist(tmp)))
          cov <- NA
      }
  } else {
      sd <- rep(NA, length(phi))
      cov <- NA
  }
  ## Output
  ans <- list(value=phi,sd=sd,cov=cov,par.fixed=par.fixed,
              cov.fixed=Vtheta,pdHess=pdHess,
              gradient.fixed=gradient.fixed)
  ## ======== Calculate bias corrected random effects estimates if requested
  if(bias.correct && !is.null(r)) {
      epsilon <- rep(0,length(phi))
      names(epsilon) <- names(phi)
      parameters <- obj$env$parameters
      parameters$TMB_epsilon_ <- epsilon ## Appends to list without changing attributes
      doEpsilonMethod <- function(chunk = NULL) {
          if(!is.null(chunk)) { ## Only do *chunk*
              mapfac <- rep(NA, length(phi))
              mapfac[chunk] <- chunk
              parameters$TMB_epsilon_ <- updateMap(parameters$TMB_epsilon_,
                                                   factor(mapfac) )
          }
          obj3 <- MakeADFun(obj$env$data,
                            parameters,
                            random = obj$env$random,
                            checkParameterOrder = FALSE,
                            DLL = obj$env$DLL,
                            silent = obj$env$silent)
          ## Get good initial parameters
          obj3$env$start <- c(par, epsilon)
          obj3$env$random.start <- expression(start[random])
          ## Test if Hessian pattern is un-changed
          h <- obj$env$spHess(random=TRUE)
          h3 <- obj3$env$spHess(random=TRUE)
          pattern.unchanged <- identical(h@i,h3@i) & identical(h@p,h3@p)
          ## If pattern un-changed we can re-use symbolic Cholesky:
          if(pattern.unchanged){
              if(!obj$env$silent)
                  cat("Re-using symbolic Cholesky\n")
              obj3$env$L.created.by.newton <- L
          } else {
              if( .Call("have_tmb_symbolic", PACKAGE = "TMB") )
                  runSymbolicAnalysis(obj3)
          }
          if(!is.null(chunk)) epsilon <- epsilon[chunk]
          par.full <- c(par.fixed, epsilon)
          i <- (1:length(par.full)) > length(par.fixed) ## epsilon indices
          grad <- obj3$gr(par.full)
          Vestimate <-
              if(bias.correct.control$sd) {
                  ## requireNamespace("numDeriv")
                  hess <- numDeriv::jacobian(obj3$gr, par.full)
                  -hess[i,i] + hess[i,!i] %*% Vtheta %*% hess[!i,i]
              } else
                  matrix(NA)
          estimate <- grad[i]
          names(estimate) <- names(epsilon)
          list(value=estimate, sd=sqrt(diag(Vestimate)), cov=Vestimate)
      }
      nsplit <- bias.correct.control$nsplit
      if(is.null(nsplit)) {
          split <- bias.correct.control$split
      } else {
          split <- split(seq_along(phi),
                         cut(seq_along(phi), nsplit))
      }
      if( is.null( split ) ){ ## Get all
          ans$unbiased <- doEpsilonMethod()
      } else {
          tmp <- lapply(split, doEpsilonMethod)
          m <- if (bias.correct.control$sd)
                   length(phi) else 1
          ans$unbiased <- list(value = rep(NA, length(phi)),
                               sd    = rep(NA, m),
                               cov   = matrix(NA, m, m))
          for(i in seq_along(split)) {
              ans$unbiased$value[ split[[i]] ] <- tmp[[i]]$value
              if (bias.correct.control$sd) {
                  ans$unbiased$sd   [ split[[i]] ] <- tmp[[i]]$sd
                  ans$unbiased$cov  [ split[[i]],
                                      split[[i]] ] <- tmp[[i]]$cov
              }
          }
      }
  }
  ## ======== bias correct but no random effects => nothing to do
  if(bias.correct && is.null(r)) {
      ans$unbiased <- list(value=ans$value, sd=ans$sd, cov=ans$cov)
      warning("'bias.correct' does nothing without random effects")
  }
  ## ======== Find marginal variances of all random effects i.e. phi(u,theta)=u
  if(!is.null(r)){
    if(is(L,"dCHMsuper")){ ## Required by inverse subset algorithm
      diag.term1 <- solveSubset(L=L, diag=TRUE)
      if(ignore.parm.uncertainty){
          diag.term2 <- 0
      } else {
          f <- obj$env$f
          w <- rep(0, length(par))
          if(!ADGradForward0Initialized) ADGradForward0Initialize()
          reverse.sweep <- function(i){
              w[i] <- 1
              f(par, order = 1, type = "ADGrad", rangeweight = w, doforward=0)[r]
          }
          nonr <- setdiff(seq_along(par), r)
          framework <- .Call("getFramework", PACKAGE=obj$env$DLL)
          if (framework != "TMBad")
              tmp <- sapply(nonr,reverse.sweep)
          else
              tmp <- f(par, order = 1, type = "ADGrad", keepx=nonr, keepy=r) ## TMBad only !!!
          if(!is.matrix(tmp)) ## Happens if length(r)==1
              tmp <- matrix(tmp, ncol=length(nonr) )
          A <- solve(hessian.random, tmp)
          diag.term2 <- rowSums((A %*% Vtheta)*A)
      }
      ans$par.random <- par[r]
      ans$diag.cov.random <- diag.term1 + diag.term2
      if(getJointPrecision){ ## Get V(u,theta)^-1
          if(length(par.fixed) == 0) {
              ans$jointPrecision <- hessian.random
          }
          else if (!ignore.parm.uncertainty) {
              G <- hessian.random %*% A
              G <- as.matrix(G) ## Avoid Matrix::cbind2('dsCMatrix','dgeMatrix')
              M1 <- cbind2(hessian.random,G)
              M2 <- cbind2(t(G), as.matrix(t(A)%*%G)+hessian.fixed )
              M <- rbind2(M1,M2)
              M <- forceSymmetric(M,uplo="L")
              dn <- c(names(par)[r],names(par[-r]))
              dimnames(M) <- list(dn,dn)
              p <- invPerm(c(r,(1:length(par))[-r]))
              ans$jointPrecision <- M[p,p]
          }
          else {
              warning("ignore.parm.uncertainty ==> No joint precision available")
          }
      }
    } else {
      warning("Could not report sd's of full randomeffect vector.")
    }
  }
  ## Copy a few selected members of the environment 'env'. In
  ## particular we need the 'skeleton' objects that allow us to put
  ## results back in same shape as original parameter list.
  ans$env <- new.env(parent = emptyenv())
  ans$env$parameters <- obj$env$parameters
  ans$env$random <- obj$env$random
  ans$env$ADreportDims <- obj2$env$ADreportDims
  class(ans) <- "sdreport"
  ans
}

##' Extract parameters, random effects and reported variables along
##' with uncertainties and optionally Chi-square statistics. Bias
##' corrected quantities are added as additional columns if available.
##'
##' @title summary tables of model parameters
##' @param object Output from \code{\link{sdreport}}
##' @param select Parameter classes to select. Can be any subset of
##' \code{"fixed"} (\eqn{\hat\theta}), \code{"random"} (\eqn{\hat u}) or
##' \code{"report"} (\eqn{\phi(\hat u,\hat\theta)}) using notation as
##' \code{\link{sdreport}}.
##' @param p.value Add column with approximate p-values
##' @param ... Not used
##' @return matrix
##' @method summary sdreport
##' @S3method summary sdreport
summary.sdreport <- function(object, select = c("all", "fixed", "random", "report"),
                             p.value=FALSE, ...)
{
  select <- match.arg(select, several.ok = TRUE)# *several* : e.g. c("fixed", "report")
  ## check if 'meth' (or "all") is among the 'select'ed ones :
  s.has <- function(meth) any(match(c(meth, "all"), select, nomatch=0L)) > 0L
  ans1 <- ans2 <- ans3 <- NULL
  if(s.has("fixed"))  ans1 <- cbind(object$par.fixed,  sqrt(diag(object$cov.fixed)))
  if(s.has("random")) ans2 <- cbind(object$par.random, sqrt(as.numeric(object$diag.cov.random)))
  if(s.has("report")) ans3 <- cbind(object$value,      object$sd)
  ans <- rbind(ans1, ans2, ans3)
  if(s.has("report")) {
      ans4 <- cbind("Est. (bias.correct)" = object$unbiased$value,
                    "Std. (bias.correct)" = object$unbiased$sd)
      if(!is.null(ans4))
          ans <- cbind(ans, rbind(NA * ans1, NA * ans2, ans4))
  }
  if(length(ans) && ncol(ans) > 0) {
    colnames(ans)[1:2] <- c("Estimate", "Std. Error")
    if(p.value) {
      ans <- cbind(ans, "z value"    = (z <- ans[,"Estimate"] / ans[,"Std. Error"]))
      ans <- cbind(ans, "Pr(>|z^2|)" = pchisq(z^2, df=1, lower.tail=FALSE))
    }
  } else
      warning("no or empty summary selected via 'select = %s'",
              deparse(select))
  ans
}

##' Print parameter estimates and give convergence diagnostic based on
##' gradient and Hessian.
##'
##' @title Print brief model summary
##' @param x Output from \code{\link{sdreport}}
##' @param ... Not used
##' @return NULL
##' @method print sdreport
##' @S3method print sdreport
print.sdreport <- function(x, ...)
{
  cat("sdreport(.) result\n")
  print(summary(x, "fixed"))
  if(!x$pdHess) {
    cat("Warning:\nHessian of fixed effects was not positive definite.\n")
  }
  cat("Maximum gradient component:", max(abs(x$gradient.fixed)),"\n")
  invisible(x)
}

##' Get estimated parameters or standard errors in the same shape as
##' the original parameter list.
##'
##' This function converts the selected column \code{what} of
##' \code{summary(x, select = c("fixed", "random"), ...)} to the same
##' format as the original parameter list (re-ordered as the template
##' parameter order). The argument \code{what} is partially matched
##' among the column names of the summary table. The actual match is
##' added as an attribute to the output.
##'
##' @title Convert estimates to original list format.
##' @param x Output from \code{\link{sdreport}}.
##' @param what Select what to convert (Estimate / Std. Error).
##' @param report Get AD reported variables rather than model parameters ?
##' @param ... Passed to \code{\link{summary.sdreport}}.
##' @return List of same shape as original parameter list.
##' @method as.list sdreport
##' @S3method as.list sdreport
##' @examples
##' \dontrun{
##' example(sdreport)
##'
##' ## Estimates as a parameter list:
##' as.list(rep, "Est")
##'
##' ## Std Errors in the same list format:
##' as.list(rep, "Std")
##'
##' ## p-values in the same list format:
##' as.list(rep, "Pr", p.value=TRUE)
##'
##' ## AD reported variables as a list:
##' as.list(rep, "Estimate", report=TRUE)
##'
##' ## Bias corrected AD reported variables as a list:
##' as.list(rep, "Est. (bias.correct)", report=TRUE)
##' }
as.list.sdreport <- function(x, what = "", report=FALSE, ...) {
    if (what == "") return (x)
    if (!report) {
        ans <- x$env$parameters
        random <- x$env$random
        par <- numeric(length(x$par.fixed) +
                       length(x$par.random))
        fixed <- rep(TRUE, length(par))
        if(length(random)>0)
            fixed[random] <- FALSE
        ## Possible choices
        opts <- colnames( summary(x, select = c("fixed", "random"), ...) )
        what <- match.arg(what, opts)
        if( any( fixed ) )
            par[ fixed ] <- summary(x, select = "fixed",  ...)[ , what]
        if( any(!fixed ) )
            par[!fixed ] <- summary(x, select = "random", ...)[ , what]
        ## Workaround utils::relist bug (?) for empty list items
        nonemp <- sapply(ans, function(x)length(x) > 0)
        nonempindex <- which(nonemp)
        skeleton <- as.relistable(ans[nonemp])
        li <- relist(par, skeleton)
        reshape <- function(x){
            if(is.null(attr(x,"map")))
                return(x)
            y <- attr(x,"shape")
            ## Handle special case where parameters are mapped to a fixed
            ## value
            if (what != "Estimate") {
                y[] <- NA
            }
            f <- attr(x,"map")
            i <- which(f >= 0)
            y[i] <- x[f[i] + 1L]
            y
        }
        for(i in seq(skeleton)){
            ans[[nonempindex[i]]][] <- as.vector(li[[i]])
        }
        for(i in seq(ans)){
            ans[[i]] <- reshape(ans[[i]])
        }
    } else { ## Reported variables
        ## Possible choices
        opts <- colnames( summary(x, select = "report", ...) )
        what <- match.arg(what, opts)
        par <- summary(x, select = "report",  ...)[ , what]
        skeleton <- lapply(x$env$ADreportDims,
                           function(dim) array(NA, dim))
        skeleton <- as.relistable(skeleton)
        ans <- relist(par, skeleton) ## Not keeping array dims !
        ans <- Map(array, ans, x$env$ADreportDims)
        class(ans) <- NULL
    }
    attr(ans, "check.passed") <- NULL
    attr(ans, "what") <- what
    ans
}
