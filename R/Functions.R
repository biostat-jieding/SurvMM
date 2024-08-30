
#==========================================================================#
# Functions with Well-Defined Structure ####
#==========================================================================#

#---------------------------------#
# The overall package description #
#---------------------------------#

#' SurvMM: Survival Mixture Model
#'
#' This package provides a framework that can estimate mixture models based on nonparametric heteroscedastic models.
#' project work
#' @importFrom methods is
#' @importFrom stats cov model.frame model.matrix model.response optim pchisq pnorm rbinom rexp runif sd uniroot
#' @docType package
#' @name SurvMM
NULL


#==========================================================================#
# Finite Mixture of Non-parametric Heteroscedastic Regression Models (with Cure Fraction) ####
#==========================================================================#

# the main function
#' @title Nonparametric Mixture Modeling of Survival Data
#'
#' @description Fit the finite mixture model that is composed of non-parametric heteroscedastic regression models.
#'
#' @aliases SurvMM.NHR.Fit
#'
#' @param yobs time to event of interest.
#' @param delta the censoring indicator, normally 1 = event of interest happens, and 0 = censoring.
#' @param X a matrix of covariates.
#' @param xGrids the grid points of fitting over the range of interested covariate.
#' @param parx.init the initial value of functional curves. The default is NULL.
#' @param numG the number of latent groups.
#' @param cure a logical value that indicates whether a cured latent group should be included.
#' @param bandwidth a list that specifies information concerning bandwidth. It contains two elements:
#'   \code{h} the candidate set of bandwidths;
#'   \code{nfold} the number of cross-validation folds.
#' @param kern the kernal function. The default is Epanechnikov kernel.
#' @param porder the order of polynomial regressions for initial values.
#' @param boots a list that specifies information concerning bandwidth. It contains two elements:
#'   \code{do} whether the bootstrap procedure should be conducted;
#'   \code{nboot} specifies the number of bootstrap sampling.
#' @param GoF a list that specifies information concerning bandwidth. It contains two elements:
#'   \code{do} whether the goodness-of-fit test statistic should be calculated;
#'   \code{siglevel} the significance level of testing.
#' @param error distribution function of error term
#' @param maxit the maximum number of iteration.
#' @param maxit.step the maximum number of iteration before likelihood increasing.
#' @param eps the error of convergence.
#' @param trace a logical to indicate whether the traces of model fitting should be printed.
#'
#' @details
#'    This is a function used to fit nonparametric mixture survival model based on the local likelihood technique.
#'
#' @examples
#' \donttest{
#'
#' #------------------------------------------------------------------------#
#' # Illustration based on a simulated dataset ####
#' #------------------------------------------------------------------------#
#'
#' # generate the used dataset (a mixture of two latent groups)
#'
#' # - true functional curves
#' parx.true <- list(
#'   px = function(x){c(0.2+0.4*sin(pi*x),0.8-0.4*sin(pi*x))},
#'   mx = function(x){c(2-0.5*sin(2*pi*x),0.5*cos(2*pi*x))},
#'   sx = function(x){c(0.4*exp(0.4*x^2),0.3*exp(-0.2*x))}
#' )
#'
#' # - generate the survival dataset and extract elements
#' set.seed(1)
#' sdata <- sdata.generate(N=500,parx=parx.true,cvalue=12)
#' yobs <- sdata$yobs
#' delta <- sdata$delta
#' X    <- as.matrix(sdata[,-c(1,2),drop=FALSE])
#'
#' # model fitting and do visualization
#'
#' # - fit the model via defined function
#' sol.FMM.NHR <- SurvMM.NHR.Fit(
#'   yobs      = yobs,
#'   delta     = delta,
#'   X         = X,
#'   xGrids    = seq(0,1,0.01),
#'   numG      = 2,
#'   bandwidth = list(
#'     h       = seq(0.03,0.21,0.03),
#'     nfold   = 5),
#'   porder    = c(0,3,0),
#'   boots     = list(do = FALSE, nboot = NULL),
#'   GoF       = list(do = FALSE, siglevel = NULL)
#' )
#'
#' # - extract fitted functional curves
#' parx <- sol.FMM.NHR$parx
#'
#' # do plot
#'
#' # - preparation
#' par(mfrow=c(1,3))
#'
#' # - plot fitted regression function
#' # __ for observed data points
#' plot(
#'   log(yobs[delta==1])~X[delta==1,],
#'   type="p",pch=16,cex=0.3,
#'   xlim=c(0,1),ylim=c(-3,3),
#'   xlab="x",ylab="Logarithm of Survival Time",
#'   main="Regression Function"
#' )
#' # __ proposed estimates
#' lines(parx$mx[,1]~xGrids,lty="dotdash",col="blue",lwd=2)
#' lines(parx$mx[,2]~xGrids,lty="dotdash",col="blue",lwd=2)
#' # __ true estimates
#' lines(sapply(xGrids,parx.true$mx)[1,]~xGrids)
#' lines(sapply(xGrids,parx.true$mx)[2,]~xGrids)
#'
#' # - plot fitted mixing proportion function
#' # __ true estimates
#' plot(
#'   sapply(xGrids,parx.true$px)[1,]~xGrids,
#'   type="l",lwd=2,
#'   xlim=c(0,1),ylim=c(0,1),
#'   xlab="x",ylab="Mixing Probability",
#'   main="Mixing Probability Function"
#' )
#' lines(sapply(xGrids,parx.true$px)[2,]~xGrids,lwd=2)
#' # __ proposed estimates
#' lines(parx$px[,1]~xGrids,lty="dotdash",col="blue",lwd=2)
#' lines(parx$px[,2]~xGrids,lty="dotdash",col="blue",lwd=2)
#'
#'
#' # - plot fitted standard error function
#' # __ true estimates
#' plot(
#'   sapply(xGrids,parx.true$sx)[1,]~xGrids,
#'   type="l",lwd=2,
#'   xlim=c(0,1),ylim=c(0.2,0.7),
#'   xlab="x",ylab="Variance of Error Term",
#'   main="Variance Function"
#' )
#' lines(sapply(xGrids,parx.true$sx)[2,]~xGrids,lwd=2)
#' # __ proposed estimates
#' lines(parx$sx[,1]~xGrids,lty="dotdash",col="blue",lwd=2)
#' lines(parx$sx[,2]~xGrids,lty="dotdash",col="blue",lwd=2)
#'
#' }
#'
#' @export SurvMM.NHR.Fit
SurvMM.NHR.Fit <- function(
    yobs,delta,X,
    xGrids = NULL,
    parx.init = NULL,
    numG = 2,
    cure = FALSE,
    bandwidth=list(
      h     = seq(0.03,0.21,by=0.03),
      nfold = 5
    ),
    kern = "Epanechnikov",
    porder = c(0,5,0),
    boots = list(do=TRUE,nboot=100),
    GoF   = list(do=TRUE,siglevel=0.05),
    error = "extreme",
    maxit = 80,
    maxit.step = 5,
    eps = 1e-5,
    trace = TRUE
){

  ### Preparations

  # define useful elements
  N     <- length(yobs)
  pX    <- ncol(X)
  if(is.null(xGrids)){xGrids <- seq(min(X),max(X),length.out=50)}
  numx <- length(xGrids)
  dist <- error.dist(error=error)

  ### find initial estimates

  if(trace==TRUE){cat("- Preparing Initial Estimates\n")}

  # if not provided before, fit a polynomial regression
  if(is.null(parx.init)){

    # fit the model
    X.poly <- do.call(cbind,lapply(1:max(porder),function(p){X^p}))
    ZXW <- lapply(1:3,function(j){if(porder[j]==0){return(NULL)}else{return(X.poly[,1:porder[j],drop=FALSE])}})
    poly.fit <- FMM.PHR.Para.fit(
      yobs=yobs,delta=delta,
      Z=ZXW[[1]],X=ZXW[[2]],W=ZXW[[3]],
      numG=numG,cure=cure,dist=dist,
      maxit=500,maxit.step=maxit.step,eps=eps
    )
    convergence.info.init <- poly.fit$convergence.info

    # record this estimator
    xGrids.poly  <- do.call(cbind,lapply(0:max(porder),function(p){xGrids^p}))
    pxtrans.init <- xGrids.poly[,1:(porder[1]+1),drop=FALSE]%*%poly.fit$par$gam
    px.init      <- exp(pxtrans.init)/apply(exp(pxtrans.init),1,sum)
    mx.init      <- xGrids.poly[,1:(porder[2]+1),drop=FALSE]%*%poly.fit$par$bet
    sx.init      <- exp(xGrids.poly[,1:(porder[3]+1),drop=FALSE]%*%poly.fit$par$the)
    parx.init    <- list(px = px.init,mx = mx.init,sx = sx.init)

  }else{
    convergence.info.init <- NULL
  }

  ### find the optimal bandwidth h

  # do EM algorithm (for different bandwidth)
  h.num <- length(bandwidth$h)
  if(h.num==1){

    # - define the optimal bandwidth directly
    bandwidth$h.optim <- bandwidth$h

  }else{

    if(trace==TRUE){cat("- Doing Cross-Validation to Obtain Optimal Bandwidth\n")}

    # - groups for Cross-Validation
    nfold.id <- sample(rep(1:bandwidth$nfold,ceiling(N/bandwidth$nfold))[1:N],N)
    CVlist <- lapply(1:bandwidth$nfold,function(ifold){(1:N)[nfold.id==ifold]})

    # - do Cross-Validation
    CVIC <- array(-Inf,dim=c(h.num,bandwidth$nfold))

    # - fit models on different h
    for(ih in 1:h.num){

      if(trace==TRUE){cat("  + CV Bandwidth",ih,"/",h.num,"\n")}

      try.ih <- try({

        for(ifold in 1:bandwidth$nfold){ # ifold <- 1; ih <- 1

          # - split the data into training and testing
          idx.test <- CVlist[[ifold]]
          yobs.train <- yobs[-idx.test]; delta.train <- delta[-idx.test]; X.train <- X[-idx.test,,drop=FALSE]
          yobs.test  <- yobs[idx.test];  delta.test  <- delta[idx.test];  X.test  <- X[idx.test,,drop=FALSE]

          # - fit the current model for a given h
          fit.ifold.ih <- FMM.NHR.EM(
            h=bandwidth$h[ih],numG=numG,cure=cure,dist=dist,parx.init=parx.init,
            yobs=yobs.train,delta=delta.train,X=X.train,
            xGrids=xGrids,kern=kern,
            maxit=maxit,maxit.step=maxit.step,eps=eps,
            trace=FALSE
          )

          # - calculate the current loss on test data
          CVIC[ih,ifold] <- FMM.NHR.Likelihood(
            yobs=yobs.test,delta=delta.test,X=X.test,
            parx=fit.ifold.ih$parx,xGrids=xGrids,
            numG=numG,cure=cure,dist=dist)

        }

      }, silent = T)

      # - obtain the optimal tuning parameter
      bandwidth$h.optim <- bandwidth$h[which.max(apply(CVIC,1,sum))]

    }

  }

  ### re-fit the model using the optimal bandwidth

  if(trace==TRUE){cat("- Fitting the Optimal Model\n")}

  # do EM algorithm (for the optimal bandwidth)
  fit.optim <- FMM.NHR.EM(
    h=bandwidth$h.optim,numG=numG,cure=cure,dist=dist,parx.init=parx.init,
    yobs=yobs,delta=delta,X=X,xGrids=xGrids,kern=kern,
    maxit=maxit,maxit.step=maxit.step,eps=eps,
    trace=trace
  )

  # extract need elements
  parx <- fit.optim$parx
  convergence.info <- fit.optim$convergence.info


  ### calculate goodness-of-fit statistic
  if(GoF$do==TRUE){

    # - calculate the goodness-of-fit statistic
    GoF$value <- FMM.NHR.GoF(
      yobs=yobs,delta=delta,X=X,parx=parx,xGrids=xGrids,
      numG=numG,cure=cure,dist=dist,h=bandwidth$h.optim,kern=kern)

  }

  ### bootstrap procedure for variance estimation
  if(boots$do==TRUE){

    if(trace==TRUE){cat("- Doing Bootstrap to Quantify Variation\n")}

    # preparations

    # == needed elements to generate survival time
    h.bootdist.stime <- bandwidth$h.optim#*N^(0.09)
    if(h.bootdist.stime==bandwidth$h.optim){
      fit.bootdist <- fit.optim
    }else{
      fit.bootdist <- FMM.NHR.EM(
        h=h.bootdist.stime,numG=numG,cure=cure,dist=dist,parx.init=parx,
        yobs=yobs,delta=delta,X=X,xGrids=xGrids,kern=kern,
        maxit=maxit,maxit.step=maxit.step,eps=eps,
        trace=FALSE
      )
    }
    if(numG==1 & cure==FALSE){
      pX.bootdist <- array(1,dim=c(N,1))
    }else{
      pXtrans.bootdist <- do.call(cbind,lapply(1:(numG-ifelse(cure,0,1)),function(g){
        lin.interpolation(x=as.vector(X),xaxis=xGrids,yaxis=log(fit.bootdist$parx$px[,g]/fit.bootdist$parx$px[,numG+ifelse(cure,1,0)]))$y}))
      exppXtrans.bootdist <- cbind(exp(pXtrans.bootdist),1)
      pX.bootdist <- exppXtrans.bootdist/apply(exppXtrans.bootdist,1,sum)
    }
    mX.bootdist <- do.call(cbind,lapply(1:numG,function(g){lin.interpolation(x=as.vector(X),xaxis=xGrids,yaxis=fit.bootdist$parx$mx[,g])$y}))
    sX.bootdist <- do.call(cbind,lapply(1:numG,function(g){exp(lin.interpolation(x=as.vector(X),xaxis=xGrids,yaxis=log(fit.bootdist$parx$sx[,g]))$y)}))
    parX.bootdist <- list(pX=pX.bootdist,mX=mX.bootdist,sX=sX.bootdist)
    # == needed elements to generate censoring time
    h.bootdist.ctime <- bandwidth$h.optim#*N^(0.09)
    WXX.bootdist <- do.call(cbind,lapply(1:N,function(i){
      kernx <- Kernel((X[i,]-as.vector(X))/h.bootdist.ctime,name=kern)/h.bootdist.ctime
      kernx/sum(kernx)
    })) # row is X, col is given x
    CTempX1.bootdist <- 1 - outer(yobs[delta==0],yobs,FUN=">=") %*% WXX.bootdist
    CTempX2.bootdist <- 1 - outer(yobs[delta==0],yobs,FUN=">" ) %*% WXX.bootdist
    CTempX.bootdist <- (CTempX1.bootdist / ifelse(CTempX2.bootdist==0,1,CTempX2.bootdist)) # row is y, col is x
    rm(WXX.bootdist);rm(CTempX1.bootdist);rm(CTempX2.bootdist)
    # == preliminary functions for generating times
    if(cure==TRUE){
      bootdist.stime <- function(t,u,px,mx,sx,dist){u-1+sum(px*c(dist$S((log(t)-mx)/sx),1))}
    }else{
      bootdist.stime <- function(t,u,px,mx,sx,dist){u-1+sum(px*dist$S((log(t)-mx)/sx))}
    }
    bootdist.ctime <- function(t,u,CTempx,yobs){u-1+prod(CTempx[yobs[delta==0]<=t])*(t<=max(yobs)) }

    # - start the bootstrap
    results.boots <- rep(list(list()),boots$nboot)
    iboot <- 1; nwrong <- 0
    while(iboot <= boots$nboot){

      if(nwrong > ceiling(boots$nboot*0.1)){stop("An error!")}
      if(trace==TRUE){if(iboot%%10==0){cat("  + Bootstrap",iboot,"/",boots$nboot,"\n")}}

      # - generate the current bootstrap sample
      if(cure==TRUE){
        stime.iboot <- sapply(1:N,function(i){
          if(rbinom(1,1,parX.bootdist$pX[i,numG+1])==0){
            val <- uniroot(bootdist.stime,c(0,1e6),u=runif(1,0,1-parX.bootdist$pX[i,numG+1]),
                           px=parX.bootdist$pX[i,],mx=parX.bootdist$mX[i,],sx=parX.bootdist$sX[i,],dist=dist)$root
          }else{val <- Inf}
          return(val)
        })
      }else{
        stime.iboot <- sapply(1:N,function(i){uniroot(
          bootdist.stime,c(0,1e6),u=runif(1,0,1),px=parX.bootdist$pX[i,],mx=parX.bootdist$mX[i,],sx=parX.bootdist$sX[i,],dist=dist
        )$root})
      }
      ctime.iboot <- sapply(1:N,function(i){uniroot(bootdist.ctime,c(0,1e6),u=runif(1,0,1),CTempx=CTempX.bootdist[,i],yobs=yobs)$root})
      yobs.iboot  <- pmin(stime.iboot,ctime.iboot)
      delta.iboot <- as.numeric(stime.iboot<=ctime.iboot)
      X.iboot     <- X

      # - fit the current model and store the results
      try.iboot <- try({
        fit.iboot <- FMM.NHR.EM(
          h=bandwidth$h.optim,numG=numG,cure=cure,dist=dist,parx.init=parx,
          yobs=yobs.iboot,delta=delta.iboot,X=X.iboot,
          xGrids=xGrids,kern=kern,
          maxit=maxit,maxit.step=maxit.step,eps=eps,
          trace=FALSE
        )
      }, silent = T)
      if( is(try.iboot, "try-error") == TRUE){nwrong<-nwrong+1;next}
      results.boots[[iboot]] <- fit.iboot$parx

      # - calculate the current goodness-of-fit statistic
      if(GoF$do==TRUE){
        results.boots[[iboot]]$GoF <- FMM.NHR.GoF(
          yobs=yobs.iboot,delta=delta.iboot,X=X.iboot,parx=fit.iboot$parx,
          xGrids=xGrids,numG=numG,cure=cure,dist=dist,h=bandwidth$h.optim,kern=kern)
      }

      # - prepare for next iteration
      iboot <- iboot + 1
    }

    # - get estimtes of standard errors
    px.boots  <- array(do.call(cbind,lapply(1:boots$nboot,function(iboot){results.boots[[iboot]]$px})),dim=c(numx,numG+ifelse(cure,1,0),boots$nboot))
    mx.boots  <- array(do.call(cbind,lapply(1:boots$nboot,function(iboot){results.boots[[iboot]]$mx})),dim=c(numx,numG,boots$nboot))
    sx.boots  <- array(do.call(cbind,lapply(1:boots$nboot,function(iboot){results.boots[[iboot]]$sx})),dim=c(numx,numG,boots$nboot))
    boots$parx.boots <- list(px=px.boots,mx=mx.boots,sx=sx.boots)
    boots$parx.sd   <- list(
      px = apply(px.boots,c(1,2),sd),
      mx = apply(mx.boots,c(1,2),sd),
      sx = apply(sx.boots,c(1,2),sd)
    )
    boots$parx.CIlow <- list(
      px = apply(px.boots,c(1,2),function(x){quantile(x,0.025)}),
      mx = apply(mx.boots,c(1,2),function(x){quantile(x,0.025)}),
      sx = apply(sx.boots,c(1,2),function(x){quantile(x,0.025)})
    )
    boots$parx.CIhig <- list(
      px = apply(px.boots,c(1,2),function(x){quantile(x,0.975)}),
      mx = apply(mx.boots,c(1,2),function(x){quantile(x,0.975)}),
      sx = apply(sx.boots,c(1,2),function(x){quantile(x,0.975)})
    )

    # - calculate critical value for goodness-of-fit statistic
    if(GoF$do==TRUE){
      GoF$value.boots <- do.call(rbind,lapply(1:boots$nboot,function(iboot){results.boots[[iboot]]$GoF}))
      GoF$critical <- apply(GoF$value.boots,2,function(x){quantile(x,1-GoF$siglevel)})
      GoF$pvalue   <- apply(t(GoF$value.boots)>=GoF$value,1,mean)
    }

  }

  # calculate the information criterion
  # - degree of freedom
  temp1 <- Kernel(0,name=kern)-0.5*integrate(f=function(u,name){Kernel(u,name)^2},lower=-Inf,upper=Inf,name=kern)$value
  temp2 <- integrate(f=Vectorize(function(u,name){
    (Kernel(u,name) - 0.5*integrate(f=(function(u2,name,u){
      Kernel(u2,name)*Kernel(u-u2,name)
    }),lower=-Inf,upper=Inf,name=name,u=u)$value)^2
  }),lower=-Inf,upper=Inf,name=kern)$value
  degree.freedom <- (3*numG-ifelse(cure,0,1))*(temp1/temp2)*diff(range(X))*temp1/bandwidth$h.optim
  # - likelihood function
  loglik <- FMM.NHR.Likelihood(yobs=yobs,delta=delta,X=X,parx=parx,xGrids=xGrids,numG=numG,cure=cure,dist=dist)
  # - final information criterion
  IC.AIC <- -2*loglik + 2*degree.freedom
  IC.BIC <- -2*loglik + log(N)*degree.freedom
  IC <- list(AIC = IC.AIC,BIC = IC.BIC)

  ### extract output values
  rownames(parx$px) <- rownames(parx$mx) <- rownames(parx$sx) <- xGrids
  colnames(parx$px) <- paste("Group",1:(numG+ifelse(cure,1,0)),sep="")
  colnames(parx$mx) <- colnames(parx$sx) <- paste("Group",1:numG,sep="")
  out <- list(
    parx = parx,
    IC = IC,
    boots = boots,
    GoF = GoF,
    bandwidth=bandwidth,
    convergence.info = c(convergence.info,list(info.init=convergence.info.init)),
    parx.init = parx.init
  )

}


# EM algorithm for fitting mixture of NHR models: fix G and h
FMM.NHR.EM <- function(
    h,numG,cure,dist,parx.init,yobs,delta,X,xGrids,kern,maxit,maxit.step,eps,trace
){

  # define useful elements
  N <- length(yobs)
  pX <- ncol(X)
  yobsln <- log(yobs)
  XI <- cbind(1,X)
  numx <- length(xGrids)

  # define old elements (needed in iterations)
  gamx.old <- gamx <- log(parx.init$px[,-(numG+ifelse(cure,1,0)),drop=FALSE]/parx.init$px[,numG+ifelse(cure,1,0)])
  betx.old <- betx <- parx.init$mx
  thex.old <- thex <- log(parx.init$sx)

  # start the iterations
  numit <- 1
  repeat{

    ## E-step: calculate the expectation of latent variable

    # - do linear interpolation
    if(numG==1 & cure==FALSE){
      pXtrans.linint <- array(0,dim=c(N,0))
    }else{
      pXtrans.linint <- do.call(cbind,lapply(1:(numG-ifelse(cure,0,1)),function(g){
        lin.interpolation(x=as.vector(X),xaxis=xGrids,yaxis=gamx.old[,g])$y
      }))
    }
    mX.linint <- do.call(cbind,lapply(1:numG,function(g){
      lin.interpolation(x=as.vector(X),xaxis=xGrids,yaxis=betx.old[,g])$y
    }))
    sX.linint <- do.call(cbind,lapply(1:numG,function(g){
      exp(lin.interpolation(x=as.vector(X),xaxis=xGrids,yaxis=thex.old[,g])$y)
    }))

    # - prepare basic elements
    exppXtrans <- cbind(exp(pXtrans.linint),1)
    mixprob <- exppXtrans/apply(exppXtrans,1,sum,na.rm=TRUE)
    residual.scaled <- (yobsln-mX.linint)/sX.linint
    hsub <- dist$h(residual.scaled)/(sX.linint*yobs)
    Ssub <- dist$S(residual.scaled)
    liksub <- (hsub^delta)*Ssub
    if(cure==TRUE){
      liksub <- cbind(liksub,1-delta)
    }

    # - the expectation of latent variables
    mixprobliksub <- mixprob*liksub
    meanlatent <- mixprobliksub/apply(mixprobliksub,1,sum,na.rm=TRUE)

    ## M-step: optimize parameters

    # - solve for each point among grids
    for( ix in 1:numx ){ # ix <- 5

      # - define the kernel weights for the current x
      Kernx.ix <- as.vector(Kernel((X-xGrids[ix])/h,name=kern)/h)

      # - update gamx for the current x (not for G = 1)
      if(numG>=2 | cure==TRUE){
        px.ix <- colSums(meanlatent*Kernx.ix)/sum(Kernx.ix,na.rm=TRUE)
        px.ix <- ifelse(px.ix==1,1-1e-6,ifelse(px.ix==0,1e-6,px.ix))
        gamx[ix,] <- log(px.ix[-(numG+ifelse(cure,1,0))]/px.ix[numG+ifelse(cure,1,0)])
      }

      # - update bet

      # = try the optimization problem
      for(g in 1:numG){

        # - fit the parametric aft model for g-th coordinate
        try.betxthex.ix.g <- try({
          numit.g <- 1
          meanlatent.g.Kernx.ix <- meanlatent[,g]*Kernx.ix
          betx.g.old <- betx.old[ix,g]; thex.g.old <- thex.old[ix,g]
          resiscal.g.old <- (yobsln-betx.g.old)/exp(thex.g.old)
          lik.g.old <- sum(meanlatent.g.Kernx.ix*(delta*dist$logh(resiscal.g.old)-delta*thex.g.old+dist$logS(resiscal.g.old)),na.rm=TRUE)
          repeat{
            # - prepare basic elements
            expthex <- exp(thex.g.old)
            resiscal <- (yobsln-betx.g.old)/expthex
            h.resiscal <- dist$h(resiscal)
            hd1.resiscal <- dist$hd1(resiscal)
            hd2.resiscal <- dist$hd2(resiscal)
            temp1 <- delta*hd1.resiscal/h.resiscal-h.resiscal
            temp2 <- delta*(hd2.resiscal*h.resiscal-hd1.resiscal^2)/(h.resiscal^2)-hd1.resiscal
            Score.betx <- -sum(meanlatent.g.Kernx.ix*temp1/expthex,na.rm=TRUE)
            Score.thex <- -sum(meanlatent.g.Kernx.ix*(temp1*resiscal+delta),na.rm=TRUE)
            InfoM.betx <- sum(meanlatent.g.Kernx.ix*temp2/(expthex^2),na.rm=TRUE)
            InfoM.thex <- sum(meanlatent.g.Kernx.ix*(temp2*(resiscal^2)+temp1*resiscal),na.rm=TRUE)
            dev.betx <- Score.betx/InfoM.betx
            dev.thex <- Score.thex/InfoM.thex
            # - update old elements
            steplen <- 1
            numit.temp <- 1
            repeat{
              betx.g.temp <- betx.g.old - steplen*dev.betx
              thex.g.temp <- thex.g.old - steplen*dev.thex
              resiscal.g.temp <- (yobsln-betx.g.temp)/exp(thex.g.temp)
              lik.g.temp <- sum(meanlatent.g.Kernx.ix*(delta*dist$logh(resiscal.g.temp)-delta*thex.g.temp+dist$logS(resiscal.g.temp)),na.rm=TRUE)
              if((numit.temp<maxit.step) & (lik.g.old>=lik.g.temp)){
                steplen <- steplen/2
                numit.temp <- numit.temp + 1
              }else{
                break
              }
            }
            # - no appropriate step, not update
            if(numit.temp>=maxit.step){
              betx.g <- betx.g.old; thex.g <- thex.g.old; lik.g <- lik.g.old
            }else{
              betx.g <- betx.g.temp; thex.g <- thex.g.temp; lik.g <- lik.g.temp
            }
            # - determine convergence
            if( max(abs(c(betx.g-betx.g.old,thex.g-thex.g.old)))>eps & numit.g<maxit ){
              betx.g.old <- betx.g; thex.g.old <- thex.g; lik.g.old <- lik.g
              numit.g <- numit.g + 1
            }else{
              break
            }
          }

        }, silent = T)

        # = record the updated result
        if(is(try.betxthex.ix.g, "try-error") == FALSE){
          betx[ix,g] <- betx.g
          thex[ix,g] <- thex.g
        }

      }

    }

    ## update or stop
    emresidual <- max(c(abs(gamx.old-gamx),abs(betx.old-betx),abs(thex.old-thex)))
    if(trace==TRUE){if(numit%%5==0){cat("  + EM Round",numit,"| Residual",round(emresidual,7),"\n")}}
    if( emresidual>eps & numit < maxit){
      gamx.old <- gamx; betx.old <- betx; thex.old <- thex
      numit <- numit + 1
    }else{
      break
    }

  }

  # transform back
  exppxtrans <- cbind(exp(gamx),1)
  px <- exppxtrans/apply(exppxtrans,1,sum,na.rm=TRUE)
  mx <- betx
  sx <- exp(thex)

  # judge the convergence
  convergence.info <- list(convergence=(numit<maxit),numit=numit,residual=emresidual)

  # collect the results and output
  out <- list(
    parx = list(px = px,mx = mx,sx = sx),
    convergence.info = convergence.info
  )
  return(out)

}

FMM.NHR.GoF <- function(yobs,delta,X,parx,xGrids,numG,cure=FALSE,dist,h,kern){

  # prepare
  N <- length(yobs)
  maxyobs <- max(yobs)

  # prepare basic elements
  if(numG==1 & cure==FALSE){pX <- array(1,dim=c(N,1))}else{
    pXtrans <- do.call(cbind,lapply(1:(numG-ifelse(cure,0,1)),function(g){lin.interpolation(x=as.vector(X),xaxis=xGrids,yaxis=log(parx$px[,g]/parx$px[,numG+ifelse(cure,1,0)]))$y}))
    exppXtrans <- cbind(exp(pXtrans),1)
    pX <- exppXtrans/apply(exppXtrans,1,sum,na.rm=TRUE)
  }
  mX  <- do.call(cbind,lapply(1:numG,function(g){lin.interpolation(x=as.vector(X),xaxis=xGrids,yaxis=parx$mx[,g])$y}))
  sX  <- do.call(cbind,lapply(1:numG,function(g){exp(lin.interpolation(x=as.vector(X),xaxis=xGrids,yaxis=log(parx$sx[,g]))$y)}))

  # calculate GoF statistic - dH
  if(cure==TRUE){
    SH0 <- sapply(1:N,function(j){mean(rowSums(pX*cbind(dist$S((log(yobs[j])-mX)/sX),1)))})
  }else{
    SH0 <- sapply(1:N,function(j){mean(rowSums(pX*dist$S((log(yobs[j])-mX)/sX)))})
  }
  SH1    <- KM.fit(tm=yobs,yobs=yobs,delta=delta)
  GoF.T1 <- sum((SH0-SH1)^2,na.rm=TRUE)

  # calculate GoF statistic - dSKM
  SH1.C <- KM.fit(tm=yobs,yobs=yobs,delta=1-delta,type="left")
  SH1.C <- ifelse(SH1.C==0,1,SH1.C)
  GoF.T2 <- sum(((SH0-SH1)^2)*delta/SH1.C,na.rm=TRUE)

  # output
  return(c(
    T1 = GoF.T1,
    T2 = GoF.T2
  ))

}


FMM.NHR.Pred.Stx <- function(tm,x,parx,xGrids,numG,cure=FALSE,error="extreme",subgroups=FALSE){

  dist <- error.dist(error=error)

  if(numG==1 & cure==FALSE){
    pXtrans.linint <- array(0,dim=c(1,0))
  }else{
    pXtrans.linint <- do.call(cbind,lapply(1:(numG-ifelse(cure,0,1)),function(g){
      lin.interpolation(x=x,xaxis=xGrids,yaxis=log(parx$px[,g]/parx$px[,numG+ifelse(cure,1,0)]))$y
    }))
  }
  exppXtrans.linint <- cbind(exp(pXtrans.linint),1)
  pX.linint <- exppXtrans.linint/apply(exppXtrans.linint,1,sum,na.rm=TRUE)

  mX.linint <- do.call(cbind,lapply(1:numG,function(g){
    lin.interpolation(x=x,xaxis=xGrids,yaxis=parx$mx[,g])$y
  }))

  sX.linint <- do.call(cbind,lapply(1:numG,function(g){
    exp(lin.interpolation(x=x,xaxis=xGrids,yaxis=log(parx$sx[,g]))$y)
  }))

  if(subgroups==FALSE){
    if(cure==TRUE){
      Stx <- sapply(tm,function(itm){
        sum(pX.linint*cbind(dist$S((log(itm)-mX.linint)/sX.linint),1),na.rm=TRUE)
      })
    }else{
      Stx <- sapply(tm,function(itm){
        sum(pX.linint*dist$S((log(itm)-mX.linint)/sX.linint),na.rm=TRUE)
      })
    }
  }else{
    Stx <- do.call(cbind,lapply(1:numG,function(g){
      sapply(tm,function(itm){dist$S((log(itm)-mX.linint[,g])/sX.linint[,g])})
    }))
  }

  return(Stx)

}


FMM.NHR.Pred.ftx <- function(tm,x,parx,xGrids,numG,cure=FALSE,error="extreme",subgroups=FALSE){

  dist <- error.dist(error=error)

  if(numG==1 & cure==FALSE){
    pXtrans.linint <- array(0,dim=c(1,0))
  }else{
    pXtrans.linint <- do.call(cbind,lapply(1:(numG-ifelse(cure,0,1)),function(g){
      lin.interpolation(x=x,xaxis=xGrids,yaxis=log(parx$px[,g]/parx$px[,numG+ifelse(cure,1,0)]))$y
    }))
  }
  exppXtrans.linint <- cbind(exp(pXtrans.linint),1)
  pX.linint <- exppXtrans.linint/apply(exppXtrans.linint,1,sum,na.rm=TRUE)

  mX.linint <- do.call(cbind,lapply(1:numG,function(g){
    lin.interpolation(x=x,xaxis=xGrids,yaxis=parx$mx[,g])$y
  }))

  sX.linint <- do.call(cbind,lapply(1:numG,function(g){
    exp(lin.interpolation(x=x,xaxis=xGrids,yaxis=log(parx$sx[,g]))$y)
  }))

  if(subgroups==FALSE){
    if(cure==TRUE){
      ftx <- sapply(tm,function(itm){
        sum(pX.linint*cbind(dist$f((log(itm)-mX.linint)/sX.linint)/(itm*sX.linint),0),na.rm=TRUE)
      })
    }else{
      ftx <- sapply(tm,function(itm){
        sum(pX.linint*dist$f((log(itm)-mX.linint)/sX.linint)/(itm*sX.linint),na.rm=TRUE)
      })
    }
  }else{
    ftx <- do.call(cbind,lapply(1:numG,function(g){
      sapply(tm,function(itm){dist$f((log(itm)-mX.linint[,g])/sX.linint[,g])/(itm*sX.linint[,g])})
    }))
  }

  return(ftx)

}



FMM.NHR.Likelihood <- function(yobs,delta,X,parx,xGrids,numG,cure=FALSE,dist){

  # - calculate basic elements
  if(numG==1 & cure==FALSE){
    pXtrans.linint <- array(0,dim=c(length(yobs),0))
  }else{
    pXtrans.linint <- do.call(cbind,lapply(1:(numG-ifelse(cure,0,1)),function(g){
      lin.interpolation(x=as.vector(X),xaxis=xGrids,yaxis=log(parx$px[,g]/parx$px[,numG+ifelse(cure,1,0)]))$y
    }))
  }
  exppXtrans.linint <- cbind(exp(pXtrans.linint),1)
  pX.linint <- exppXtrans.linint/apply(exppXtrans.linint,1,sum,na.rm=TRUE)

  mX.linint <- do.call(cbind,lapply(1:numG,function(g){
    lin.interpolation(x=as.vector(X),xaxis=xGrids,yaxis=parx$mx[,g])$y
  }))

  sX.linint <- do.call(cbind,lapply(1:numG,function(g){
    exp(lin.interpolation(x=as.vector(X),xaxis=xGrids,yaxis=log(parx$sx[,g]))$y)
  }))

  residual.scaled <- (log(yobs)-mX.linint)/sX.linint
  hsub   <- dist$h(residual.scaled)/(sX.linint*yobs)
  Ssub   <- dist$S(residual.scaled)
  liksub <- (hsub^delta)*Ssub
  if(cure==TRUE){liksub <- cbind(liksub,1-delta)}
  loglik <- sum(log(apply(pX.linint*liksub,1,sum)),na.rm=TRUE)

  # collect the results and output
  return(loglik)

}



#=================================================================================#
# A General Function for Fitting Finite Mixture of Parametric AFT Models  ####
#=================================================================================#

# finite mixture of parametric aft models
FMM.PHR.Para.fit <- function(
    yobs,delta,X,Z=NULL,W=NULL,
    cure=FALSE,
    numG=NULL,
    dist=NULL,
    maxit=300,maxit.step=5,eps=1e-5,
    params.init=NULL
){

  ### Preparations

  # basic definitions
  N <- length(yobs)
  pX <- ncol(X)
  pZ <- ifelse(is.null(Z),0,ncol(Z))
  pW <- ifelse(is.null(W),0,ncol(W))
  if(is.null(colnames(X))){colnames(X)<-paste("X",1:pX,sep="")}
  if(is.null(colnames(Z))){if(pZ>0){colnames(Z)<-paste("Z",1:pZ,sep="")}}
  if(is.null(colnames(W))){if(pW>0){colnames(W)<-paste("W",1:pW,sep="")}}
  yobsln <- log(yobs)
  XI <- as.matrix(cbind(rep(1,N),X))
  ZI <- as.matrix(cbind(rep(1,N),Z))
  WI <- as.matrix(cbind(rep(1,N),W))

  ### do EM algorithm

  # initial estimates
  if(is.null(params.init)){
    gam.old <- gam <- array(0,dim=c(pZ+1,numG-ifelse(cure,0,1)))
    bet.old <- bet <- array(0,dim=c(pX+1,numG))
    the.old <- the <- array(0,dim=c(pW+1,numG))
    the.old[1,] <- the[1,] <- seq(-0.3,0.3,length.out=numG)
    # gam.old[1,] <- gam[1,] <- seq(-0.2,0.2,length.out=numG-ifelse(cure,0,1))
  }else{
    if(cure){
      gam.old <- gam <- params.init$gam
    }else{
      gam.old <- gam <- params.init$gam[,-numG,drop=FALSE]
    }
    gam.old <- gam <- params.init$gam[,-numG,drop=FALSE]
    bet.old <- bet <- params.init$bet
    the.old <- the <- params.init$the
  }

  # start the iterations
  numit <- 1
  repeat{

    ## E-step: calculate the expectation of latent variable

    # - calculate basic quantities
    expZIgam <- cbind(exp(ZI%*%gam.old),1)
    mixprob  <- expZIgam/apply(expZIgam,1,sum,na.rm=TRUE)
    expWIthe <- exp(WI%*%the.old)
    residual.scaled <- (yobsln-XI%*%bet.old)/expWIthe
    hsub   <- dist$h(residual.scaled)/(expWIthe*yobs)
    Ssub   <- dist$S(residual.scaled)
    liksub <- (hsub^delta)*Ssub
    if(cure==TRUE){liksub <- cbind(liksub,1-delta)}

    # - the expectation of latent variables
    mixprobliksub <- mixprob*liksub
    meanlatent <- mixprobliksub/apply(mixprobliksub,1,sum,na.rm=TRUE)


    ## M-step: optimize parameters

    # - update gam
    if(numG>=2 | cure==TRUE){

      # = try the optimization problem
      try.gam <- try({

        numit.outer <- 1
        gam.temp <- gam.old
        lik.old <- sum(meanlatent*log(expZIgam/apply(expZIgam,1,sum,na.rm=TRUE)),na.rm=TRUE)
        repeat{
          # - iteratively solve each coordinate
          for(g in 1:(numG-ifelse(cure,0,1))){
            # - solve the optimization problem for g-th coordinate
            numit.inner.g <- 1
            repeat{
              # - calculate basic elements
              expZIgam <- cbind(exp(ZI%*%gam),1)
              mixprob.g  <- expZIgam[,g]/apply(expZIgam,1,sum,na.rm=TRUE)
              Score.g <- apply(ZI*(meanlatent[,g] - mixprob.g),2,sum,na.rm=TRUE)
              InfoM.g <- - t(ZI*mixprob.g*(1-mixprob.g))%*%ZI
              dev <- as.vector(MASS::ginv(InfoM.g)%*%Score.g)
              # - update old elements
              steplen <- 1
              numit.temp <- 1
              repeat{
                gam.inner.temp <- gam
                gam.inner.temp[,g] <- gam[,g] - steplen*dev
                expZIgam.temp <- cbind(exp(ZI%*%gam.inner.temp),1)
                lik.temp <- sum(meanlatent*log(expZIgam.temp/apply(expZIgam.temp,1,sum,na.rm=TRUE)),na.rm=TRUE)
                if((numit.temp<maxit.step) & (lik.old>=lik.temp)){
                  steplen <- steplen/2
                  numit.temp <- numit.temp + 1
                }else{
                  break
                }
              }
              # - no appropriate step, not update
              if(numit.temp>=maxit.step){
                gam.inner.g <- gam[,g]; lik <- lik.old
              }else{
                gam.inner.g <- gam.inner.temp[,g]; lik <- lik.temp
              }
              # - judge the convergence for this coordinate
              if( max(abs(dev))>eps & numit.inner.g<maxit ){
                gam[,g] <- gam.inner.g; lik.old <- lik
                numit.inner.g <- numit.inner.g + 1
              }else{
                gam[,g] <- gam.inner.g; lik.old <- lik
                break
              }
            }
          }
          # - judge the convergence for this coordinate
          if( max(abs(gam.temp-gam))>eps & numit.outer<maxit ){
            gam.temp <- gam
            numit.outer <- numit.outer + 1
          }else{
            break
          }
        }

      }, silent = T)

      # = record the updated result
      if(is(try.gam,"try-error")==TRUE){
        gam <- gam.old
      }

    }


    # - update bet
    for(g in 1:numG){

      # - fit the parametric aft model for g-th coordinate

      # = try the optimization problem
      try.betthe.g <- try({

        numit.g <- 1
        meanlatent.g <- meanlatent[,g]
        bet.g.old <- bet.old[,g]; the.g.old <- the.old[,g]
        WIthe.g.old <- as.vector(WI%*%the.g.old)
        resiscal.g.old <- as.vector((yobsln-XI%*%bet.g.old)/exp(WIthe.g.old))
        lik.g.old <- sum(meanlatent.g*(delta*dist$logh(resiscal.g.old)-delta*WIthe.g.old+dist$logS(resiscal.g.old)),na.rm=TRUE)
        repeat{
          # - prepare basic elements
          expWIthe <- as.vector(exp(WI%*%the.g.old))
          resiscal <- as.vector((yobsln-XI%*%bet.g.old)/expWIthe)
          h.resiscal <- dist$h(resiscal)
          hd1.resiscal <- dist$hd1(resiscal)
          hd2.resiscal <- dist$hd2(resiscal)
          temp1 <- delta*hd1.resiscal/h.resiscal-h.resiscal
          temp2 <- delta*(hd2.resiscal*h.resiscal-hd1.resiscal^2)/(h.resiscal^2)-hd1.resiscal
          Score.bet <- -apply(XI*meanlatent.g*temp1/expWIthe,2,sum,na.rm=TRUE)
          Score.the <- -apply(WI*meanlatent.g*(temp1*resiscal+delta),2,sum,na.rm=TRUE)
          InfoM.bet <- t(XI*meanlatent.g*temp2/(expWIthe^2))%*%XI
          InfoM.the <- t(WI*meanlatent.g*(temp2*(resiscal^2)+temp1*resiscal))%*%WI
          dev.bet <- as.vector(MASS::ginv(InfoM.bet)%*%Score.bet)
          dev.the <- as.vector(MASS::ginv(InfoM.the)%*%Score.the)
          # - update old elements
          steplen <- 1
          numit.temp <- 1
          repeat{
            bet.g.temp <- bet.g.old - steplen*dev.bet
            the.g.temp <- the.g.old - steplen*dev.the
            WIthe.temp <- as.vector(WI%*%the.g.temp)
            resiscal.g.temp <- as.vector((yobsln-XI%*%bet.g.temp)/exp(WIthe.temp))
            lik.g.temp <- sum(meanlatent.g*(delta*dist$logh(resiscal.g.temp)-delta*WIthe.temp+dist$logS(resiscal.g.temp)),na.rm=TRUE)
            if((numit.temp<maxit.step) & (lik.g.old>=lik.g.temp)){
              steplen <- steplen/2
              numit.temp <- numit.temp + 1
            }else{
              break
            }
          }
          # - no appropriate step, not update
          if(numit.temp>=maxit.step){
            bet.g <- bet.g.old; the.g <- the.g.old; lik.g <- lik.g.old
          }else{
            bet.g <- bet.g.temp; the.g <- the.g.temp; lik.g <- lik.g.temp
          }
          # - determine convergence
          if( max(abs(c(bet.g-bet.g.old,the.g-the.g.old)))>eps & numit.g<maxit ){
            bet.g.old <- bet.g; the.g.old <- the.g; lik.g.old <- lik.g
            numit.g <- numit.g + 1
          }else{
            break
          }
        }

      }, silent = T)

      # = record the updated result
      if(is(try.betthe.g,"try-error")==FALSE){
        bet[,g] <- bet.g
        the[,g] <- the.g
      }

    }

    ## update or stop
    emresidual <- max(c(abs(bet.old-bet),abs(gam.old-gam),abs(the.old-the)))
    if( emresidual>eps & numit<maxit ){
      bet.old <- bet; gam.old <- gam; the.old <- the
      numit <- numit + 1
    }else{
      break
    }

  }

  # judge the convergence
  convergence.info <- list(convergence=(numit<maxit),numit=numit,residual=emresidual)

  ### tidy the results
  nameG <- paste("Group",1:(numG+ifelse(cure,1,0)),sep="")
  nameX <- c('Intercept',colnames(X))
  nameZ <- c('Intercept',colnames(Z))
  nameW <- c('Intercept',colnames(W))
  par <- list(
    gam = array(cbind(gam,0),dim=c(pZ+1,numG+ifelse(cure,1,0)),dimnames=list(nameZ,nameG)),
    bet = array(bet,dim=c(pX+1,numG),dimnames=list(nameX,nameG[1:numG])),
    the = array(the,dim=c(pW+1,numG),dimnames=list(nameW,nameG[1:numG]))
  )

  ### extract output values
  out <- list(
    par = par,
    convergence.info = convergence.info
  )

}


#=================================================================================#
# Other Useful Functions  ####
#=================================================================================#

# specify basic functions of a given distribution
error.dist <- function(error){

  # extreme, gaussian, logistic
  if(error=="gaussian"){
    # - the error term follows the standard normal distribution
    dist <- list(
      S    = function(t){return(1-pnorm(t))},
      f    = function(t){return(dnorm(t))},
      h    = function(t){return(dnorm(t)/(1-pnorm(t)))},
      hd1  = function(t){h=dnorm(t)/(1-pnorm(t));return(-t*h+h^2)},
      hd2  = function(t){h=dnorm(t)/(1-pnorm(t));return((t^2-1)*h-3*t*h^2+2*h^3)},
      logS = function(t){return(log(1-pnorm(t)))},
      logf = function(t){return(log(dnorm(t)))},
      logh = function(t){return(log(dnorm(t)/(1-pnorm(t))))},
    )
  }else if(error=="extreme"){
    # - the error term follows the standard extreme value distribution
    dist <- list(
      S    = function(t){exp(-exp(t))},
      f    = function(t){exp(t-exp(t))},
      h    = function(t){exp(t)},
      hd1  = function(t){exp(t)},
      hd2  = function(t){exp(t)},
      logS = function(t){-exp(t)},
      logf = function(t){t-exp(t)},
      logh = function(t){t}
    )
  }else if(error=="logistic"){
    # - the error term follows the standard logistic distribution
    dist <- list(
      S    = function(t){return(1/(1+exp(t)))},
      f    = function(t){S=1/(1+exp(t));return(S-S^2)},
      h    = function(t){return(1-1/(1+exp(t)))},
      hd1  = function(t){S=1/(1+exp(t));return(S-S^2)},
      hd2  = function(t){S=1/(1+exp(t));return(-S+3*S^2-2*S^3)},
      logS = function(t){return(-log(1+exp(t)))},
      logf = function(t){S=1/(1+exp(t));return(log(S-S^2))},
      logh = function(t){return(log(1-1/(1+exp(t))))}
    )
  }

  return(dist)

}

# define different kernel functions
Kernel <- function(u, name="Epanechnikov"){
  # the kernel function (not re-scaled) with only one continuous variable

  if(name=="Uniform"){
    ret <- 0.5*(abs(u)<=1)
  }else if(name=="Gaussian"){
    ret <- exp(-u^2/2)/sqrt(2*pi)
  }else if(name=="Epanechnikov"){
    ret <- 0.75*(1-u^2)*(abs(u)<=1)
  }else if(name=="Quartic"){
    ret <- (15/16)*(1-u^2)^2*(abs(u)<=1)
  }
  return(ret)
} # Kernel(0.2,name="Gaussian")

# do linear interpolation
lin.interpolation <- function(x,xaxis,yaxis){

  # preparations
  nx <- length(x)
  naxis <- length(xaxis)

  # do linear interpolation
  y <- sapply(x,function(xc){
    idx <- sum(xaxis <= xc)+1
    if(idx<=2){
      val <- yaxis[1]+(xc-xaxis[1])*(yaxis[2]-yaxis[1])/(xaxis[2]-xaxis[1])
    }else if(idx>=naxis){
      val <- yaxis[naxis-1]+(xc-xaxis[naxis-1])*(yaxis[naxis]-yaxis[naxis-1])/(xaxis[naxis]-xaxis[naxis-1])
    }else{
      val <- yaxis[idx-1]+(xc-xaxis[idx-1])*(yaxis[idx]-yaxis[idx-1])/(xaxis[idx]-xaxis[idx-1])
    }
    val
  })

  # output
  return(list(
    x=x,y=y
  ))

}

# calculate Kaplan-Meier estimator
KM.fit <- function(tm,yobs,delta,type="right"){

  ### preparation ###
  N <- length(yobs)
  y.order <- order(yobs)
  y.sort <- yobs[y.order]
  delta.sort <- delta[y.order]
  yobs.1max <- max(yobs[delta==1])

  # calculate the values of KM at pre-specified time points tm
  prods <- 1-delta.sort/(N-(1:N)+1)
  if(type=="right"){
    KMt <- sapply(tm,function(tmi){prod(prods[y.sort<=tmi])}) # right-continuous
  }else{
    KMt <- sapply(tm,function(tmi){prod(prods[y.sort<tmi])}) # left-continuous
  }

  # output
  return(KMt)

}

#=================================================================================#
# Data Generating Function  ####
#=================================================================================#
sdata.generate <-  function(N,parx,cvalue,cure=FALSE){

  # generate covariates
  X <- array(NA, dim=c(N,1))
  X[,1] <- runif(N,0,1)

  # generating survival time
  stime <- rep(NA, N)
  for(i in 1:N){
    pxi <- parx$px(X[i,])
    gid <- which(rmultinom(1,size=1,prob=pxi)==1)
    if(cure==TRUE & gid==length(pxi)){
      stime[i] <- Inf
    }else{
      mxi <- parx$mx(X[i,])[gid]
      sxi <- parx$sx(X[i,])[gid]
      stime[i] <- exp(mxi + sxi*log(rexp(1)))
    }

  }

  # generate censoring time
  ctime <- rexp(N, 1/cvalue)
  delta <- as.numeric(stime<=ctime)
  censor.rate <- 1-mean(delta)
  cure.rate   <- mean(stime==Inf)

  # delta and observed failure time
  yobs <- pmin(stime,ctime)

  # info
  info <- list(
    censor.rate = censor.rate,
    cure.rate   = cure.rate
  ); info
  print(info)

  # output
  out <- data.frame(yobs=yobs,delta=delta,X)
  return(out)

}






