##' Fits a binary mixture of linear models in which the probability of
##' class membership is related to the covariates through a probit
##' regression model.
##'
##' The model assumes the observations are drawn from a two component
##' mixture, where each component is described by a different linear
##' model.  The probability that an individual observation is a member
##' of one component or the other is modelled by a probit regression.
##'
##' The model is fit by Gibbs sampling, assuming uniform priors for
##' the regression coefficients of the two linear model and the probit
##' regression, and a (common) Gamma prior for the precision (inverse
##' variance) of the two linear models.  The probit component of the
##' model is sampled by the method of Albert and Chib.
##'
##' @title Binary Mixture of Linear Models
##' @param formula1 An object of class formula: a symbolic description
##'   of the linear model for the first component
##' @param formula2 An object of class formula: a symbolic description
##'   of the linear model for the second component
##' @param formulap An object of class formula: a symbolic description
##'   of the probit model for the probability an observation is
##' @param data A dataframe containing the variables from the model.
##' @param nsamp Number of samples to draw.
##' @param nthin The thinning rate.
##' @param tau.prior Parameters of the (common) Gamma prior for the
##'   precision of the two models.
##' @param start A list of initial values for sigma, betap and b
##' @return An object of class \code{bmixlm} with columns
##' \item{call}{the matched call}
##' \item{nsamp}{the number of samples retained after thinning}
##' \item{beta1}{matrix of samples of the coefficients of the first
##' linear model}
##' \item{beta2}{matrix of samples of the coefficients of the second
##' linear model}
##' \item{betap}{matrix of samples of the coefficients of the probit
##' model}
##' \item{sigma}{two column matrix of samples of the standard
##' deviations of the errors for the two models}
##' \item{data}{the input dataframe}
##' \item{pW}{effective degrees of freedom for the fitted model}
##' \item{WAIC}{WAIC for the fitted model}
##' \item{restart}{final sigma, betap and b for restart purposes}
##' @references
##'   Albert, J. H., & Chib, S. (1993). Bayesian analysis of
##'   binary and polychotomous response data. Journal of the American
##'   statistical Association, 88(422), 669-679.
##' @importFrom stats dnorm pnorm qnorm rnorm rbinom rgamma runif model.frame model.response model.matrix
##' @export
bmixlm <- function(formula1,formula2,formulap,data,
                   nsamp=1000,nthin=3,
                   tau.prior=c(0.01,0.01),
                   start=list(sigma=c(1.0E-4,1.0E-4))) {

  cl <- match.call()

  ## Truncated Normal deviates
  rztnorm <- function(n,mu,y) {
    cs <- pnorm(0,mu)
    qnorm(runif(n,ifelse(y==0,0,cs),ifelse(y==0,cs,1)),mu)
  }

  ## Gibbs sampler for the linear model
  gibbs.reg <- function(y,X,sigma) {
    V <- solve(crossprod(X))
    beta <- drop(V%*%crossprod(X,y) + t(chol(V))%*%rnorm(ncol(X),0,sigma))
    r <- y-X%*%beta
    tau <- rgamma(1,tau.prior[1]+length(r)/2,tau.prior[2]+crossprod(r)/2)
    list(beta=beta,sigma=1/sqrt(tau))
  }

  ## Extract responses and model matrices for the three model
  ## components.
  mf <- model.frame(formula1,data)
  y1 <- model.response(mf)
  X1 <- model.matrix(formula1,mf)

  mf <- model.frame(formula2,data)
  y2 <- model.response(mf)
  X2 <- model.matrix(formula2,mf)

  mf <- model.frame(formulap,data)
  Xp <- model.matrix(formulap,mf)

  if(any(y1!=y2) || is.null(y1)) stop("Component formulae should have identical non-null responses")

  ## Initialize storage
  Beta1 <- matrix(0,nsamp,ncol(X1),dimnames=list(NULL,colnames(X1)))
  Beta2 <- matrix(0,nsamp,ncol(X2),dimnames=list(NULL,colnames(X2)))
  Betap <- matrix(0,nsamp,ncol(Xp),dimnames=list(NULL,colnames(Xp)))
  Sigma <- matrix(0,nsamp,2,dimnames=list(NULL,c("sigma1","sigma2")))
  B <- matrix(0,nrow(data),nsamp)
  Q <- matrix(0,nrow(data),nsamp)
  L <- matrix(0,nrow(data),nsamp)

  ## Initialize
  n <- nrow(data)
  Pp <- solve(crossprod(Xp))
  Lp <- t(chol(Pp))
  Pp <- Pp%*%t(Xp)

  sigma <- if(is.null(start$sigma)) c(1.0E-4,1.0E-4) else start$sigma
  betap <- if(is.null(start$betap)) rep_len(0,ncol(Xp)) else start$betap
  etap <- Xp%*%betap
  ## If b not supplied, restart from expected value of q
  b <- if(is.null(start$b)) rbinom(n,1,pnorm(etap)) else start$b

  ## Gibbs sample
  for(k1 in seq_len(nsamp)) {
    for(k2 in seq_len(nthin)) {

      ## Probit regression (Albert and Chib)
      z <- rztnorm(n,etap,b)
      betap <- drop(Pp%*%z+Lp%*%rnorm(ncol(Xp)))
      etap <- Xp%*%betap
      p <- pnorm(etap)

      ## First linear model
      if(sum(b==0)>ncol(X1)) {
        fit <- gibbs.reg(y1[b==0],X1[b==0,,drop=FALSE],sigma[1])
        beta1 <- fit$beta
        sigma[1] <- fit$sigma
      }

      ## Second linear model
      if(sum(b==1)>ncol(X2)) {
        fit <- gibbs.reg(y2[b==1],X2[b==1,,drop=FALSE],sigma[2])
        beta2 <- fit$beta
        sigma[2] <- fit$sigma
      }

      ## Mixing fractions
      f1 <- (1-p)*dnorm(y1,X1%*%beta1,sigma[1])
      f2 <- p*dnorm(y2,X2%*%beta2,sigma[2])
      q <- f2/(f1+f2)
      b <- rbinom(n,1,q)

    }
    ## Record samples
    Beta1[k1,] <- beta1
    Beta2[k1,] <- beta2
    Betap[k1,] <- betap
    Sigma[k1,] <- sigma
    B[,k1] <- b
    L[,k1] <- f1+f2
  }

  ## WAIC calculations
  lpd <- sum(log(rowMeans(L)))
  L <- log(L)
  pW <- sum(apply(L,1,var))

  structure(
    list(call=cl,nsamp=nsamp,
         beta1=Beta1,
         beta2=Beta2,
         betap=Betap,
         sigma=Sigma,
         b=B,
         pW=pW,
         WAIC=-2*(lpd-pW),
         data=data,
         restart=list(sigma=sigma,betap=betap,b=b)),
    class="bmixlm")
}

##' Print method for objects of class \code{bmixlm}
##'
##' @title Printing bmixlm fits
##' @param x An object of class \code{\link{bmixlm}}
##' @param ... Currently ignored
##' @export
print.bmixlm <- function(x,...) {
  cl <- x$call
  if(!is.null(cl$start)) cl$start <- NULL
  cat("\nCall:\n",paste(deparse(cl),sep="\n",collapse="\n"),
      "\nSamples: ",x$nsamp,"\n\n",sep ="")
}

##' Computes a list of summary statistics for the \code{bmixlm} model
##' fit \code{object}.
##'
##' @title Summarizing \code{bmixlm} Fits
##' @param object An object of class \code{\link{bmixlm}}
##' @param x An object of class \code{\link{bmixlm}}
##' @param digits The number of significant digits to use when printing.
##' @param ... Currently ignored
##' @return Returns an object of class \code{summary.bmixlm}, with components
##' \item{\code{call}}{The original \code{bmixlm} call}
##' \item{\code{nsamp}}{The number of Gibbs samples drawn}
##' \item{\code{beta1}}{Summary table for the coefficients for the linear model for the first component}
##' \item{\code{beta2}}{Summary table for the coefficients for the linear model for the second component}
##' \item{\code{betap}}{Summary table for the coefficients for the probit model}
##' \item{\code{sigma}}{Summary table for the standard deviations for the two linear models}
##' \item{\code{vcov1}}{Variance covariance matrix of the coefficients for the linear model for the first component}
##' \item{\code{vcov2}}{Variance covariance matrix of the coefficients for the linear model for the second component}
##' \item{\code{vcovp}}{Variance covariance matrix of the coefficients for the probit model}
##' \item{\code{pW}}{WAIC effective number of parameters}
##' \item{\code{WAIC}}{WAIC}
##' @importFrom stats sd quantile var
##' @export
summary.bmixlm <- function(object,...) {
  smry <- function(ps)
    t(apply(ps,2,function(p) c(`Mean`=mean(p),`Std Dev`=sd(p),quantile(p,c(0.5,0.025,0.975)))))

  structure(
    list(
      call=object$call,
      nsamp=object$nsamp,
      beta1=smry(object$beta1),
      beta2=smry(object$beta2),
      betap=smry(object$betap),
      sigma=smry(object$sigma),
      vcov1=var(object$beta1),
      vcov2=var(object$beta2),
      vcovp=var(object$betap),
      pW=object$pW,
      WAIC=object$WAIC),
    class="summary.bmixlm")
}

##' @rdname summary.bmixlm
##' @importFrom stats printCoefmat
##' @export
print.summary.bmixlm <- function(x,digits=max(3L,getOption("digits")-3L),...) {
  cl <- x$call
  if(!is.null(cl$start)) cl$start <- NULL
  cat("\nCall:\n",paste(deparse(cl),sep="\n",collapse="\n"),
      "\nSamples: ",x$nsamp,"\n\n",sep ="")
  cat("\nComponent 1:\n")
  printCoefmat(x$beta1, digits=digits)
  cat("\nComponent 2:\n")
  printCoefmat(x$beta2, digits=digits)
  cat("\nProbit:\n")
  printCoefmat(x$betap, digits=digits)
  cat("\nErrors:\n")
  printCoefmat(x$sigma, digits=digits)
  cat("\nWAIC:\n")

  print(data.frame(pW = x$pW, WAIC = x$WAIC,row.names = "WAIC"), ...)
}



##' Updates and optionally refits a \code{\link{bmixlm}} fit.
##'
##' If the \code{betap} and \code{sigma} arguments are not specified,
##' they are determined from the existing fit.
##'
##' @title Update and Refit a \code{bmixlm} Model
##' @param object An object of class \code{\link{bmixlm}}
##' @param ... Arguments to update
##' @param evaluate If true evaluate the new call else return the call.
##' @return If \code{evaluate = TRUE} the fitted object, otherwise the updated call.
##' @importFrom stats getCall update
##' @export
update.bmixlm <- function (object,...,evaluate=TRUE) {
    call <- getCall(object)
    extras <- match.call(expand.dots=FALSE)$...
    restart <- object$restart
    if(!is.na(match("formulap",names(extras)))) {
      restart$betap <- NULL
      restart$b <- NULL
    }
    if(is.na(match("start",names(extras)))) extras$start <- restart
    if(length(extras)) {
      existing <- !is.na(match(names(extras), names(call)))
      for(a in names(extras)[existing])
        call[[a]] <- extras[[a]]
      if(any(!existing))
        call <- as.call(c(as.list(call), extras[!existing]))
    }
    if(evaluate) eval(call, parent.frame()) else call
}


##' Extract a set of model coefficients from a fitted \code{bmixlm} object.
##'
##' The fitted object contains of four sets of parameters: the
##' coefficients for the two component models, the coefficients for
##' the probit model, and the error standard deviations from the two
##' component models. The \code{which} argument determines which of
##' these parameter sets is returned:
##' \describe{
##'   \item{\code{"comp1"}}{coefficients for the first component model}
##'   \item{\code{"comp2"}}{coefficients for the second component model}
##'   \item{\code{"probit"}}{coefficients for the probit binary model}
##'   \item{\code{"error"}}{error standard deviations for the two components}
##' }
##'
##' @title Coefficients of a \code{bmixlm} Object
##' @param object An object of class \code{bmixlm}
##' @param which The parameter set to extract (see details).
##' @param type Whether to return the posterior mean or samples from the posterior.
##' @param ... Currently unused.
##' @return If \code{type="mean"} the coefficients are returned as a
##'   vector, and if \code{type="samples"} the coefficients are
##'   returned as an array of samples from the posterior
##' @export
coef.bmixlm <- function(object,
                        which=c("probit","comp1","comp2","error"),
                        type=c("mean","samples"),...) {
  which <- match.arg(which)
  type <- match.arg(type)
  cf <- switch(which,
               probit=object$betap,
               comp1=object$beta1,
               comp2=object$beta2,
               error=object$sigma)
  if(type=="mean") colMeans(cf) else cf
}

##' Extract a model matrix for one component of a \code{bmixlm} object.
##'
##' A \code{bmixlm} object depends upon three model matrices: the
##' matrices for the two component models and the model matrix for
##' the probit model. The \code{which} argument determines which of
##' these matrices is returned:
##' \describe{
##'   \item{\code{"comp1"}}{the model matrix for the first component model}
##'   \item{\code{"comp2"}}{the model matrix for the second component model}
##'   \item{\code{"probit"}}{the model matrix for the probit binary model}
##' }
##'
##' If a dataframe is supplied via the \code{data} argument, it is
##' used to construct the model matrix, otherwise the model matrix is
##' constructed from the data used to generate the fitted object
##'
##' @title Model matrices for a \code{bmixlm} Object
##' @param object An object of class \code{bmixlm}
##' @param which The parameter set to extract (see details).
##' @param data A data frame or \code{NULL}.
##' @param ... Currently unused.
##' @return The model matrix for the selected model component.
##' @export
model.matrix.bmixlm <- function(object,which=c("probit","comp1","comp2"),data=NULL,...) {
  which <- match.arg(which)
  if(is.null(data)) data <- object$data
  formula <- switch(which,
                    probit=as.formula(object$call$formulap),
                    comp1=as.formula(object$call$formula1),
                    comp2=as.formula(object$call$formula2))
  mf <- model.frame(formula,data)
  model.matrix(formula,mf)
}




##' Return the full range of predicted values that can be computed
##' from the model when the responses are available for the prediction
##' data set.
##'
##' This function returns the predicted values that are conditional
##' only on the responses for the data to which the model is fitted:
##' \itemize{
##'   \item the predicted values for the two component linear models
##'   \item the posterior predictive probability p of membership of
##'   the second component
##' }
##' together with the predicted values conditional on the responses
##' from both the data to which the model is fitted, and the
##' prediction data set:
##' \itemize{
##'   \item the residuals or prediction error given the observed
##'   response in the prediction data
##'   \item the posterior probability q of membership of the second
##'   component,
##' }
##' The posterior predictive probability p predicts the probability
##' that a new response will be drawn from the second component, the
##' posterior probability q predicts the probability that the observed
##' response was drawn from the second component.
##'
##' If \code{type="mean"} the function returns posterior means as a
##' dataframe, otherwise it returns samples from the posterior as a
##' list of arrays.
##'
##' @title Predicted Values for a \code{bmixlm} Object
##' @param object An object of class \code{bmixlm}
##' @param type Whether to return the posterior means or samples from
##' the posterior.
##' @param data A dataframe for which predictions are required
##' @param standardize Whether to standardize the residuals.
##' @param ... Currently unused.
##' @return Returns a list or dataframe with elements
##' \item{\code{y}}{vector of responses from data}
##' \item{\code{y1}}{predicted values for the first component model}
##' \item{\code{y2}}{predicted values for the second component model}
##' \item{\code{r1}}{residuals for the first component model}
##' \item{\code{r2}}{residuals for the second component model}
##' \item{\code{p}}{posterior predictive probabilities of membership
##' of the second component}
##' \item{\code{q}}{posterior probabilities of membership of the
##' second component}
##' If \code{type="mean"} return a dataframe of posterior means is
##' returned, and if \code{type="samples"} return a list of arrays of
##' samples from the posterior.
##' @importFrom stats as.formula model.frame model.matrix model.response dnorm pnorm
##' @export
predictAll <- function(object,type=c("mean","samples"),data=NULL,standardize=FALSE,...) {
  type <- match.arg(type)

  if(is.null(data)) data <- object$data
  sigma <- object$sigma

  formula <- as.formula(object$call$formula1)
  mf <- model.frame(formula,data)
  X <- model.matrix(formula,mf)
  yhat1 <- X%*%t(object$beta1)
  y1 <- model.response(mf)
  if(!is.null(y1)) {
    r1 <- y1-yhat1
    if(standardize) r1 <- r1/rep.int(sigma[,1],rep.int(length(r1),nrow(sigma)))
  }

  formula <- as.formula(object$call$formula2)
  mf <- model.frame(formula,data)
  X <- model.matrix(formula,mf)
  yhat2 <- X%*%t(object$beta2)
  y2 <- model.response(mf)
  if(!is.null(y2)) {
    r2 <- y2-yhat2
    if(standardize) r2 <- r2/rep.int(sigma[,2],rep.int(length(r2),nrow(sigma)))
  }

  formula <- as.formula(object$call$formulap)
  mf <- model.frame(formula,data)
  X <- model.matrix(formula,mf)
  p <- pnorm(X%*%t(object$betap))
  if(!is.null(y1) && !is.null(y2)) {
    f1 <- (1-p)*dnorm(y1,yhat1,rep.int(sigma[,1],rep.int(nrow(p),nrow(sigma))))
    f2 <- p*dnorm(y2,yhat2,rep.int(sigma[,2],rep.int(nrow(p),nrow(sigma))))
    q <- f2/(f1+f2)
  }

  switch(type,
         mean=if(is.null(y1) || is.null(y2))
                data.frame(y1=rowMeans(yhat1),y2=rowMeans(yhat2),p=rowMeans(p))
              else
                data.frame(y=y1,y1=rowMeans(yhat1),y2=rowMeans(yhat2),
                           r1=rowMeans(r1),r2=rowMeans(r2),p=rowMeans(p),q=rowMeans(q)),
         samples=if(is.null(y1) || is.null(y2))
                   list(y1=yhat1,y2=yhat2,p=p)
                 else
                   list(y=y1,y1=yhat1,y2=yhat2,r1=r1,r2=r2,p=p,q=q))
}


##' Calculate the fitted values from the two component models, and the
##' probabilities of membership of the second component.
##'
##' If \code{type="mean"} the function returns posterior mean
##' quantities as a dataframe, otherwise it returns samples from the
##' posterior as a list of arrays.
##'
##' @title Fitted Values for a \code{bmixlm} Object
##' @param object An object of class \code{bmixlm}
##' @param type Whether to return the posterior mean or samples from
##'   the posterior.
##' @param ... Currently unused.
##' @return Returns a list or dataframe with elements
##' \item{\code{y1}}{predicted values for the first component model}
##' \item{\code{y2}}{predicted values for the second component model}
##' \item{\code{p}}{posterior predictive probabilities of membership
##' of the second component}
##' If \code{type="mean"} return a dataframe of posterior means is
##' returned, and if \code{type="samples"} return a list of arrays of
##' samples from the posterior.
##' @importFrom stats as.formula model.frame model.matrix pnorm
##' @export
fitted.bmixlm <- function(object,type=c("mean","samples"),...) {
  type <- match.arg(type)

  data <- object$data

  formula <- as.formula(object$call$formula1)
  mf <- model.frame(formula,data)
  X <- model.matrix(formula,mf)
  yhat1 <- X%*%t(object$beta1)

  formula <- as.formula(object$call$formula2)
  mf <- model.frame(formula,data)
  X <- model.matrix(formula,mf)
  yhat2 <- X%*%t(object$beta2)

  formula <- as.formula(object$call$formulap)
  mf <- model.frame(formula,data)
  X <- model.matrix(formula,mf)
  p <- pnorm(X%*%t(object$betap))
  b <- object$b

  switch(type,
         mean=data.frame(y1=rowMeans(yhat1),
                         y2=rowMeans(yhat2),
                         p=rowMeans(p),
                         b=rowMeans(b)),
         samples=list(y1=yhat1,
                      y2=yhat2,
                      p=p,
                      b=b))

}

##' Calculate the residuals from the two component models, and the
##' probabilities of membership of the second component.
##'
##' If \code{type="mean"} the function returns posterior mean
##' quantities as a dataframe, otherwise it returns samples from the
##' posterior as a list of arrays.
##'
##' @title Residuals for a \code{bmixlm} Object
##' @param object An object of class \code{bmixlm}
##' @param type Whether to return the posterior mean or samples from
##' the posterior.
##' @param standardize Whether to standardize the residuals.
##' @param ... Currently unused.
##' @return Returns a list or dataframe with elements
##' \item{\code{y1}}{residuals for the first component model}
##' \item{\code{y2}}{residuals for the second component model}
##' \item{\code{p}}{posterior predictive probabilities of membership
##' of the second component}
##' \item{\code{q}}{posterior probabilities of membership of the
##' second component}
##' \item{\code{b}}{binary indicators of membership of the second component,
##' conditional on the observed response}
##' If \code{type="mean"} return a dataframe of posterior means is
##' returned, and if \code{type="samples"} return a list of arrays of
##' samples from the posterior.
##' @importFrom stats as.formula model.frame model.response model.matrix pnorm
##' @export
residuals.bmixlm <- function(object,type=c("mean","samples"),standardize=FALSE,...) {
  type <- match.arg(type)

  data <- object$data
  sigma <- object$sigma

  formula <- as.formula(object$call$formula1)
  mf <- model.frame(formula,data)
  X <- model.matrix(formula,mf)
  yhat1 <- X%*%t(object$beta1)
  y1 <- model.response(mf)
  r1 <- y1-yhat1
  if(standardize) r1 <- r1/rep.int(sigma[,1],rep.int(length(r1),nrow(sigma)))

  formula <- as.formula(object$call$formula2)
  mf <- model.frame(formula,data)
  X <- model.matrix(formula,mf)
  yhat2 <- X%*%t(object$beta2)
  y2 <- model.response(mf)
  r2 <- y2-yhat2
  if(standardize) r2 <- r2/rep.int(sigma[,2],rep.int(length(r2),nrow(sigma)))

  formula <- as.formula(object$call$formulap)
  mf <- model.frame(formula,data)
  X <- model.matrix(formula,mf)
  p <- pnorm(X%*%t(object$betap))
  b <- object$b

  switch(type,
         mean=data.frame(y1=rowMeans(r1),
                         y2=rowMeans(r2),
                         p=rowMeans(p),
                         b=rowMeans(b)),
         samples=list(y1=r1,
                      y2=r2,
                      p=p,
                      b=b))
}


##' Calculate the predicted values that can be computed from a
##' fitted model.
##'
##' This function returns the predicted values conditional on the
##' original observed data and the prediction covariates only. An
##' expanded range of predictions that include quantities conditional
##' on the observed responses in the prediction data can be calculated
##' with \code{\link{predictAll}}.
##'
##' If \code{type="mean"} the function returns posterior mean
##' quantities as a dataframe, otherwise it returns samples from the
##' posterior as a list of arrays.
##'
##' @title Predicted Values for a \code{bmixlm} Object
##' @param object An object of class \code{bmixlm}
##' @param newdata A dataframe for which predictions are required
##' @param type Whether to return the posterior mean or samples
##' from the posterior.
##' @param ... Currently unused.
##' @return Returns a list or dataframe with elements
##' \item{\code{y1}}{predicted values for the first component model}
##' \item{\code{y2}}{predicted values for the second component model}
##' \item{\code{p}}{posterior predictive probabilities of membership
##' of the second component}
##' If \code{type="mean"} return a dataframe of posterior means is
##' returned, and if \code{type="samples"} return a list of arrays of
##' samples from the posterior.
##' @importFrom stats as.formula model.frame model.matrix pnorm
##' @export
predict.bmixlm <- function(object,newdata=NULL,type=c("mean","samples"),...) {
  type <- match.arg(type)

  if(is.null(newdata)) newdata <- object$data

  formula <- as.formula(object$call$formula1)
  mf <- model.frame(formula,newdata)
  X <- model.matrix(formula,mf)
  yhat1 <- X%*%t(object$beta1)

  formula <- as.formula(object$call$formula2)
  mf <- model.frame(formula,newdata)
  X <- model.matrix(formula,mf)
  yhat2 <- X%*%t(object$beta2)

  formula <- as.formula(object$call$formulap)
  mf <- model.frame(formula,newdata)
  X <- model.matrix(formula,mf)
  p <- pnorm(X%*%t(object$betap))

  switch(type,
         mean=data.frame(y1=rowMeans(yhat1),
                         y2=rowMeans(yhat2),
                         p=rowMeans(p)),
         samples=list(y1=yhat1,
                      y2=yhat2,
                      p=p))
}


##' Plot a trace plot of a set of parameters from a \code{bmixlm} object.
##'
##' The fitted object contains of four sets of parameters: the
##' coefficients for the two component models, the coefficients for
##' the probit model, and the error standard deviations from the two
##' component models. The \code{which} argument determines which of
##' these parameter sets is plotted:
##' \describe{
##'   \item{\code{"comp1"}}{coefficients for the first component model}
##'   \item{\code{"comp2"}}{coefficients for the second component model}
##'   \item{\code{"probit"}}{coefficients for the probit binary model}
##'   \item{\code{"error"}}{error standard deviations for the two components}
##' }
##'
##' @title Trace plots for \code{bmixlm} objects
##' @param x An object of class \code{bmixlm}
##' @param which The coefficient set to plot
##' @param main The main title for the plot
##' @param ... Additional options to \code{plot.ts}.
##' @importFrom graphics plot
##' @importFrom stats as.ts coef
##' @export
plot.bmixlm <- function(x,which=c("probit","comp1","comp2","error"),main=which,...) {
  which <- match.arg(which)
  plot(as.ts(coef(x,which=which,type="samples")),main=main,...)
}

##' Returns all the samples of the parameters as one large matrix.
##'
##' @title Convert \code{bmixlm} object to a matrix
##' @param x An object of class \code{bmixlm}.
##' @param ... Currently ignored.
##' @return A matrix of samples of model parameters.
##' @importFrom stats setNames
##' @export
as.matrix.bmixlm <- function(x,...) {
  r <- cbind(x$beta1,x$beta2,x$betap,x$sigma)
  colnames(r) <- c(paste("comp1",colnames(x$beta1),sep="."),
                   paste("comp2",colnames(x$beta2),sep="."),
                   paste("probit",colnames(x$betap),sep="."),
                   paste("error",colnames(x$sigma),sep="."))
  r
}


##' Predict probabilities of component membership for a new data set.
##'
##' For each row in the dataframe \code{data}, predict
##' \itemize{
##'   \item \code{p}: the posterior predictive probability of
##'   membership of the second component (conditional only on the
##'   observed covariates in the prediction data set)
##'
##'   \item \code{p}: the posterior probability of membership of the
##'   second component (conditional on both the observed response and
##'   covariates in the prediction data set)
##' }
##' If \code{type="mean"} the function returns posterior means as a
##' dataframe, otherwise it returns samples from the posterior as a
##' list of arrays.
##'
##' @title Predicted Probabilities of Component Membership
##' @param object An object of class \code{bmixlm}
##' @param type Whether to return the posterior means or samples
##' from the posterior.
##' @param data A dataframe for which predictions are required.
##' @return Returns a list or dataframe with elements
##' \item{\code{p}}{posterior predictive probabilities of membership
##' of the second component}
##' \item{\code{q}}{posterior probabilities of membership of the
##' second component}
##' If \code{type="mean"} return a dataframe of posterior means is
##' returned, and if \code{type="samples"} return a list of arrays of
##' samples from the posterior.
##' @importFrom stats as.formula model.frame model.matrix model.response dnorm pnorm
##' @export
classify <- function(object,type=c("mean","samples"),data=NULL) {
  type <- match.arg(type)

  if(is.null(data)) data <- object$data
  sigma <- object$sigma

  formula <- as.formula(object$call$formula1)
  mf <- model.frame(formula,data)
  X <- model.matrix(formula,mf)
  yhat1 <- X%*%t(object$beta1)
  y1 <- model.response(mf)

  formula <- as.formula(object$call$formula2)
  mf <- model.frame(formula,data)
  X <- model.matrix(formula,mf)
  yhat2 <- X%*%t(object$beta2)
  y2 <- model.response(mf)

  formula <- as.formula(object$call$formulap)
  mf <- model.frame(formula,data)
  X <- model.matrix(formula,mf)
  p <- pnorm(X%*%t(object$betap))
  f1 <- (1-p)*dnorm(y1,yhat1,rep.int(sigma[,1],rep.int(nrow(p),nrow(sigma))))
  f2 <- p*dnorm(y2,yhat2,rep.int(sigma[,2],rep.int(nrow(p),nrow(sigma))))
  q <- f2/(f1+f2)

  switch(type,
         mean=data.frame(p=rowMeans(p),
                         q=rowMeans(q)),
         samples=list(p=p,
                      q=q))
}

##' New observations are simulated from the posterior predictive
##' distribution.
##'
##' If \code{nsim} is not \code{NULL}, \code{update} is called to
##' generate \code{nsim} new samples from the posterior, and new
##' observations are generated from these. Otherwise new observations
##' are generated from the samples contained in the fitted object.
##'
##' @title Simulate Responses from a \code{bmixlm} Object
##' @param object An object of class \code{bmixlm}.
##' @param nsim Number of response vectors to simulate.
##' @param seed An object specifying if and how the random number
##'   generator should be initialized.
##' @param ... Additional arguments to be passed to \code{update}.
##' @return A dataframe of simulated responses.
##' @importFrom stats as.formula model.frame model.matrix runif rnorm simulate
##' @export
simulate.bmixlm <- function(object,nsim=1,seed=NULL,...) {

  if (!exists(".Random.seed",envir=.GlobalEnv,inherits=FALSE)) runif(1)
  if(is.null(seed))
    RNGstate <- get(".Random.seed",envir=.GlobalEnv)
  else {
    R.seed <- get(".Random.seed",envir=.GlobalEnv)
    set.seed(seed)
    RNGstate <- structure(seed,kind=as.list(RNGkind()))
    on.exit(assign(".Random.seed",R.seed,envir=.GlobalEnv))
  }

  object <- if(is.null(nsim)) object else update(object,nsamp=nsim,...)
  data <- object$data
  sigma <- object$sigma

  n <- nrow(data)
  m <- nrow(sigma)

  formula <- as.formula(object$call$formula1)
  mf <- model.frame(formula,data)
  X <- model.matrix(formula,mf)
  y1 <- rnorm(n*m,X%*%t(object$beta1),rep.int(sigma[,1],rep.int(n,m)))

  formula <- as.formula(object$call$formula2)
  mf <- model.frame(formula,data)
  X <- model.matrix(formula,mf)
  y2 <- rnorm(n*m,X%*%t(object$beta2),rep.int(sigma[,2],rep.int(n,m)))

  y <- (1-object$b)*y1+object$b*y2
  colnames(y) <- paste("sim",1:m,sep="_")
  as.data.frame(y)
}




