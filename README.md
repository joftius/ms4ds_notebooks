# README 

These are notebooks to accompany the course [Modern Statistics (and Causal Inference) for Data Science](http://ms4ds.com) 

We'll make use of these packages:

```
# Fairly standard packages
install.packages(c("tidyverse", "devtools", "glmnet"))
# High-dimensional / ML inference specialty packages
install.packages(c("hdi", "RPtests", "SLOPE", "stabs",
                   "knockoff", "selectiveInference", "hdm"))
devtools::install_github("swager/crossEstimation")
```