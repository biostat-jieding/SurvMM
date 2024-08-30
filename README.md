# *SurvMM* R Package

## Nonparametric mixture modeling of survival data with latent heterogeneity

This is an R package for analyzing heterogeneous **Surv**ival data based on finite **M**ixture **M**odels.
- The underlying mixture model is composed of *nonparametric heteroscedastic regression models*.

## Package description and included main functions

Installation of this package can be done locally after downloading the package manually from this github website. We will also upload this package to the Comprehensive R Archive Network (CRAN) so that it can be downloaded as a standard R package. Currently, it can be loaded using R command
```R
devtools::install_github("biostat-jieding/SurvMM")
library(SurvMM)
```

The main function included in our R package is *SurvMM.NHR.Fit()*. It can be called via the following R synopsis:
```R
SurvMM.NHR.Fit(
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
)
```
We refer to its help page for more detailed explanations of the corresponding arguments (typing *?SurvMM.NHR.Fit()*). 

## Illustration via a simulated dataset

An example of the use of this package can be found in the following of this file.

### Data preparation

Set the true underlying functional curves:
```R
parx.true <- list( 
  px = function(x){c(0.2+0.4*sin(pi*x),0.8-0.4*sin(pi*x))},
  mx = function(x){c(2-0.5*sin(2*pi*x),0.5*cos(2*pi*x))},
  sx = function(x){c(0.4*exp(0.4*x^2),0.3*exp(-0.2*x))} 
)
```

Generate the survival dataset (function *sdata.generate()* can be found within our R package):
```R
set.seed(1)
sdata <- sdata.generate(N=500,parx=parx.true,cvalue=12)
yobs <- sdata$yobs
delta <- sdata$delta
X    <- as.matrix(sdata[,-c(1,2),drop=FALSE])
```

### Model fitting 

Fit the model via provided function:
```R
xGrids <- seq(0,1,0.01)
sol.FMM.NHR <- SurvMM.NHR.Fit(
  yobs      = yobs, 
  delta     = delta, 
  X         = X,
  xGrids    = xGrids,
  numG      = 2,
  bandwidth = list(
    h       = seq(0.03,0.21,0.03),
    nfold   = 5),
  porder    = c(0,3,0),
  boots     = list(do = FALSE, nboot = NULL),
  GoF       = list(do = FALSE, siglevel = NULL)
)
```

Extract fitted functional curves:
```R
parx <- sol.FMM.NHR$parx
```

### Visualization

Preparation：
```R
par(mfrow=c(1,3))
```

Plot fitted regression function:
```R
# __ for observed data points
plot(
  log(yobs[delta==1])~X[delta==1,],
  type="p",pch=16,cex=0.3,
  xlim=c(0,1),ylim=c(-3,3),
  xlab="x",ylab="Logarithm of Survival Time",
  main="Regression Function"
)
# __ proposed estimates
lines(parx$mx[,1]~xGrids,lty="dotdash",col="blue",lwd=2)
lines(parx$mx[,2]~xGrids,lty="dotdash",col="blue",lwd=2)
# __ true estimates
lines(sapply(xGrids,parx.true$mx)[1,]~xGrids)
lines(sapply(xGrids,parx.true$mx)[2,]~xGrids)
```

Plot fitted mixing proportion function:
```R
# __ true estimates
plot(
  sapply(xGrids,parx.true$px)[1,]~xGrids,
  type="l",lwd=2,
  xlim=c(0,1),ylim=c(0,1),
  xlab="x",ylab="Mixing Probability",
  main="Mixing Probability Function"
)
lines(sapply(xGrids,parx.true$px)[2,]~xGrids,lwd=2)
# __ proposed estimates
lines(parx$px[,1]~xGrids,lty="dotdash",col="blue",lwd=2)
lines(parx$px[,2]~xGrids,lty="dotdash",col="blue",lwd=2)
```

Plot fitted standard error function:
```R
# __ true estimates
plot(
  sapply(xGrids,parx.true$sx)[1,]~xGrids,
  type="l",lwd=2,
  xlim=c(0,1),ylim=c(0.2,0.7),
  xlab="x",ylab="Variance of Error Term",
  main="Variance Function"
)
lines(sapply(xGrids,parx.true$sx)[2,]~xGrids,lwd=2)
# __ proposed estimates
lines(parx$sx[,1]~xGrids,lty="dotdash",col="blue",lwd=2)
lines(parx$sx[,2]~xGrids,lty="dotdash",col="blue",lwd=2)
```

*This R package was contributed by **Jie Ding** and **Ingrid Van Keilegom**.*
