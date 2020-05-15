causalML
================
Joshua Loftus
5/3/2020

``` r
library(tidyverse)
```

    ## ── Attaching packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────────────────── tidyverse 1.3.0 ──

    ## ✓ ggplot2 3.3.0     ✓ purrr   0.3.4
    ## ✓ tibble  3.0.1     ✓ dplyr   0.8.5
    ## ✓ tidyr   1.0.2     ✓ stringr 1.4.0
    ## ✓ readr   1.3.1     ✓ forcats 0.5.0

    ## ── Conflicts ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────── tidyverse_conflicts() ──
    ## x dplyr::filter() masks stats::filter()
    ## x dplyr::lag()    masks stats::lag()

``` r
library(glmnet)
```

    ## Loading required package: Matrix

    ## 
    ## Attaching package: 'Matrix'

    ## The following objects are masked from 'package:tidyr':
    ## 
    ##     expand, pack, unpack

    ## Loaded glmnet 4.0

``` r
library(broom)
```

Two stage least squares
-----------------------

2SLS: consider the model

*Y* = *D**θ* + *X**β* + *ϵ*<sub>*Y*</sub>
 and assume that *D* has not been randomized, but that *X* are observed confounders satisfying

*D* = *X**α* + *ϵ*<sub>*D*</sub>
 and the errors are strictly exogenous.

By linearity, we could also fit

*Y* = *X*(*θ**α* + *β*)+(*θ**ϵ*<sub>*D*</sub> + *ϵ*<sub>*Y*</sub>)=*X**γ* + *ϵ*<sub>*Y*</sub>′
 but our goal is to estimate *θ*.

-   Regress *Y* on *D* and *X*
-   Regress *Y* on *X*, then regress the residuals on *D*
-   Regress *Y* on *X*, regress *D* on *X*, then regress the *Y* residuals on the *D* residuals
-   Regress *D* on *X*, then regress the *Y* on the predicted $\\hat D$
-   Regress *D* on *X*, then regress the *Y* on the *D* residuals

``` r
instance <- function(n = 50, x0 = 0, sx = 1, d0 = 0, alpha = 0, sd = 1, y0 = 0, theta = 0, beta = 0, sy = 1) {
  X <- x0 + sqrt(sx) * rnorm(n)
  D <- d0 + alpha * X + sqrt(sd) * rnorm(n)
  Y <- theta * D + beta * X + sqrt(sy) * rnorm(n)  
  
  lmYDX <- lm(Y ~ D + X)
  lmDX <- lm(D~X)
  lmYX <- lm(Y~X)
  Yres <- resid(lmYX)
  Ypred <- predict(lmYX)
  Ypres <- Y - coef(lmYDX)[3] * X
  Dres <- resid(lmDX)
  Dpred <- predict(lmDX)
  
  output <- c(
    coef(lmYDX)[2],
    coef(lm(Ypres ~ D))[2],
    coef(lm(Yres ~ Dres))[2],
    coef(lm(Y ~ Dres))[2],
    coef(lm(Y ~ Dpred))[2], 
    coef(lm(Ypred ~ Dpred))[2],
    cov(Y,X)/cov(D,X),
    coef(lm(Y ~ D))[2],
    coef(lm(Yres ~ D))[2]
  )
  names(output) <- c(
    "Y~D+X", "Ypres~D", "Yres~Dres", "Y~Dres", "Y~Dhat", "Yhat~Dhat", "Wald", "Y~D(OVB)", "Yres(OVB)~D"
  )
  output
}

MCMSE <- function(nsim = 1000, n = 50, x0 = 0, sx = 1, d0 = 0, alpha = 0, sd = 1, y0 = 0, theta = 0, beta = 0, sy = 1) {
  rowMeans(replicate(nsim, instance(n, x0, sx, d0, alpha, sd, y0, theta, beta, sy))) - theta
}
```

``` r
instance(n = 50, alpha = 2, theta = 5, beta = 3, sd = 1)
```

    ##       Y~D+X     Ypres~D   Yres~Dres      Y~Dres      Y~Dhat   Yhat~Dhat 
    ##    5.143370    5.143370    5.143370    5.143370    6.226994    6.226994 
    ##        Wald    Y~D(OVB) Yres(OVB)~D 
    ##    6.226994    6.013256    1.014495

``` r
instance(n = 500, alpha = 2, theta = 5, beta = 3, sd = 1)
```

    ##       Y~D+X     Ypres~D   Yres~Dres      Y~Dres      Y~Dhat   Yhat~Dhat 
    ##    5.050432    5.050432    5.050432    5.050432    6.499634    6.499634 
    ##        Wald    Y~D(OVB) Yres(OVB)~D 
    ##    6.499634    6.211757    1.003244

``` r
instance(n = 5000, alpha = 2, theta = 5, beta = 3, sd = 1)
```

    ##       Y~D+X     Ypres~D   Yres~Dres      Y~Dres      Y~Dhat   Yhat~Dhat 
    ##    5.002372    5.002372    5.002372    5.002372    6.487787    6.487787 
    ##        Wald    Y~D(OVB) Yres(OVB)~D 
    ##    6.487787    6.182793    1.027116

``` r
MCMSE(nsim = 1000, n = 50, alpha = 2, theta = 5, beta = 3, sd = 1)
```

    ##       Y~D+X     Ypres~D   Yres~Dres      Y~Dres      Y~Dhat   Yhat~Dhat 
    ## -0.00102508 -0.00102508 -0.00102508 -0.00102508  1.50175867  1.50175867 
    ##        Wald    Y~D(OVB) Yres(OVB)~D 
    ##  1.50175867  1.19744153 -3.99749919

``` r
MCMSE(nsim = 1000, n = 100, alpha = 2, theta = 5, beta = 3, sd = 1)
```

    ##         Y~D+X       Ypres~D     Yres~Dres        Y~Dres        Y~Dhat 
    ##  0.0008659043  0.0008659043  0.0008659043  0.0008659043  1.5048676859 
    ##     Yhat~Dhat          Wald      Y~D(OVB)   Yres(OVB)~D 
    ##  1.5048676859  1.5048676859  1.2021372040 -3.9975598799

``` r
MCMSE(nsim = 1000, n = 500, alpha = 2, theta = 5, beta = 3, sd = 1)
```

    ##        Y~D+X      Ypres~D    Yres~Dres       Y~Dres       Y~Dhat    Yhat~Dhat 
    ##  0.001020753  0.001020753  0.001020753  0.001020753  1.500123474  1.500123474 
    ##         Wald     Y~D(OVB)  Yres(OVB)~D 
    ##  1.500123474  1.200000556 -3.999595164

``` r
# IV case: beta = 0
instance(n = 100, alpha = 2, theta = 5, beta = 0, sd = 1)
```

    ##       Y~D+X     Ypres~D   Yres~Dres      Y~Dres      Y~Dhat   Yhat~Dhat 
    ##   4.9214473   4.9214473   4.9214473   4.9214473   5.0534325   5.0534325 
    ##        Wald    Y~D(OVB) Yres(OVB)~D 
    ##   5.0534325   5.0295032   0.8922716

``` r
#MCMSE(nsim = 1000, n = 100, alpha = 2, theta = 5, beta = 0, sd = 1)
```

Non-linear case
---------------

Higher-dimensional case
-----------------------

``` r
correlated_gaussian_design <- function(n, p, rho) {
  x <- matrix(rnorm(n*p), nrow = n)
  if (rho == 0) return(x)
  z <- matrix(rep(t(rnorm(n)), p), nrow = n)
  sqrt(1-rho)*x + sqrt(rho)*z
}

cvlasso <- function(x, y) {
  cvfit <- cv.glmnet(x, y, intercept = FALSE)
  lfit <- glmnet(x, y, intercept = FALSE)
  minlambda <- cvfit$lambda.min
  coef(lfit, s = minlambda)[-1]
}

hd_instance <- function(n = 100, p = 200, rho = 0, alpha_sparsity = 0, alpha = 0, theta = 0, beta_sparsity = 0, beta = 0) {
  X <- correlated_gaussian_design(n, p, rho)
  D <- rnorm(n)
  al <- rep(0, p)
  if (alpha_sparsity > 0) {
    al[1:alpha_sparsity] <- alpha
    D <- D + X %*% al
  }
  Y <- theta * D + rnorm(n)  
  be <- rep(0, p)
  if (beta_sparsity > 0) {
    be[1:beta_sparsity] <- beta
    Y <- Y + X %*% be
  }
  
  DX <- cbind(D, X)
  lassoYDX <- cvlasso(DX, Y)
  lassoDX <- cvlasso(X, D)
  lassoYX <- cvlasso(X, Y)
  
  Ypred <- X %*% lassoYX
  Yres <- Y - Ypred
  Ypres <- Y - X %*% lassoYDX[-1]
  Dpred <- X %*% lassoDX
  Dres <- D - Dpred
  
  # DML
  nsplits <- 2
  splits <- matrix(sample(1:n), nrow = nsplits)
  theta_hats <- rep(0, nsplits)
  DMLtheta_hats <- rep(0, nsplits)
  for (j in 1:nsplits) {
    D_test <- splits[j,]
    D_train <- setdiff(1:n, D_test)
    stage1ML <- cvlasso(DX[D_train,], Y[D_train])
    DML_Ypres <- Y - X %*% stage1ML[-1]
    
    theta_hats[j] <- mean(D[D_test] * DML_Ypres[D_test])/mean(D[D_test]^2)
      #coef(lm(DML_Ypres[D_test] ~ D[D_test]))[2]
    
    lassoDX_train <- cvlasso(X[D_train, ], D[D_train])
    Dpred_train <- X %*% lassoDX_train
    Dres_train <- D - Dpred_train
    
    DMLtheta_hats[j] <- mean(Dres_train[D_test] * DML_Ypres[D_test])/mean(Dres_train[D_test] * D[D_test])
  }
  
  output <- c(
    mean(theta_hats),
    mean(DMLtheta_hats),
    lassoYDX[1],
    coef(lm(Ypres ~ D))[2],
    coef(lm(Yres ~ Dres))[2],
    coef(lm(Y ~ Dres))[2],
    coef(lm(Y ~ Dpred))[2], 
    coef(lm(Ypred ~ Dpred))[2],
    coef(lm(Y ~ D))[2],
    coef(lm(Yres ~ D))[2]
  )
  names(output) <- c(
    "2SML", "DML", "Y~D+X", "Ypres~D", "Yres~Dres", "Y~Dres", "Y~Dhat", "Yhat~Dhat", "Y~D(OVB)", "Yres(OVB)~D"
  )
  output
}

hd_MCMSE <- function(nsim = 1000, n = 100, p = 200, rho = 0, alpha_sparsity = 0, alpha = 0, theta = 0, beta_sparsity = 0, beta = 0) {
  rowMeans( (replicate(nsim, hd_instance(n, p, rho, alpha_sparsity, alpha, theta, beta_sparsity, beta) - theta)^2 )) 
}
```

``` r
# does it work?
hd_instance(alpha_sparsity = 5, alpha = 1, theta = 3.14, beta_sparsity = 10, beta = 1)
```

    ##        2SML         DML       Y~D+X     Ypres~D   Yres~Dres      Y~Dres 
    ##    3.991128    3.729357    3.832809    3.934258    3.116553    7.656147 
    ##      Y~Dhat   Yhat~Dhat    Y~D(OVB) Yres(OVB)~D 
    ##    5.004080    4.411949    4.171102    0.757446

``` r
hd_instance(n = 400, alpha_sparsity = 5, alpha = 1, theta = 3.14, beta_sparsity = 10, beta = 1)
```

    ##        2SML         DML       Y~D+X     Ypres~D   Yres~Dres      Y~Dres 
    ##    3.500587    3.289996    3.290774    3.313372    2.897293    4.951408 
    ##      Y~Dhat   Yhat~Dhat    Y~D(OVB) Yres(OVB)~D 
    ##    4.419945    4.181978    3.886989    0.628835

``` r
hd_instance(n = 800, alpha_sparsity = 5, alpha = 1, theta = 3.14, beta_sparsity = 10, beta = 1)
```

    ##        2SML         DML       Y~D+X     Ypres~D   Yres~Dres      Y~Dres 
    ##   3.3022724   3.1753201   3.2252271   3.2414375   3.0835798   4.2544419 
    ##      Y~Dhat   Yhat~Dhat    Y~D(OVB) Yres(OVB)~D 
    ##   4.3459497   4.1445312   3.9356644   0.6263294

``` r
#hd_MCMSE(alpha_sparsity = 5, alpha = 1, theta = 3.14, beta_sparsity = 10, beta = 1)
#hd_MCMSE(n = 400, alpha_sparsity = 5, alpha = 1, theta = 3.14, beta_sparsity = 10, beta = 1)
#hd_MCMSE(n = 800, alpha_sparsity = 5, alpha = 1, theta = 3.14, beta_sparsity = 10, beta = 1)
```
