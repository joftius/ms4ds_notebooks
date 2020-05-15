lassoinf
================
Joshua Loftus
4/16/2020

Submodels vs full models
------------------------

``` r
n <- 100
Z <- rnorm(n)
X1 <- rnorm(n)
X2 <- 2*X1 + Z
Y <- 3*X1 + 1.5*X2 + rnorm(n)

summary(lm(Y~X1+X2))
```

    ## 
    ## Call:
    ## lm(formula = Y ~ X1 + X2)
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -2.2445 -0.6507  0.0307  0.6315  2.0827 
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept) -0.03067    0.09240  -0.332    0.741    
    ## X1           3.37041    0.23215  14.518   <2e-16 ***
    ## X2           1.27917    0.11022  11.606   <2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.8888 on 97 degrees of freedom
    ## Multiple R-squared:  0.9802, Adjusted R-squared:  0.9797 
    ## F-statistic:  2395 on 2 and 97 DF,  p-value: < 2.2e-16

``` r
summary(lm(Y~X1))
```

    ## 
    ## Call:
    ## lm(formula = Y ~ X1)
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -3.4854 -0.7193  0.1867  0.9288  3.0919 
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)  -0.2830     0.1381  -2.049   0.0431 *  
    ## X1            5.8726     0.1323  44.375   <2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 1.367 on 98 degrees of freedom
    ## Multiple R-squared:  0.9526, Adjusted R-squared:  0.9521 
    ## F-statistic:  1969 on 1 and 98 DF,  p-value: < 2.2e-16

``` r
c(cor(Y,X1), cor(Y,X2))
```

    ## [1] 0.9760076 0.9679985

``` r
summary(lm(Y~X2))
```

    ## 
    ## Call:
    ## lm(formula = Y ~ X2)
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -4.1112 -0.8786  0.0668  0.9431  4.5375 
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)  0.19320    0.16146   1.197    0.234    
    ## X2           2.76532    0.07242  38.185   <2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 1.575 on 98 degrees of freedom
    ## Multiple R-squared:  0.937,  Adjusted R-squared:  0.9364 
    ## F-statistic:  1458 on 1 and 98 DF,  p-value: < 2.2e-16
