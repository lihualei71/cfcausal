# cfcausal

## Overview
This package implements the weighted conformal inference procedure for counterfactuals and individual treatment effects proposed in the paper: [Conformal Inference of Counterfactuals and Individual Treatment Effects](https://arxiv.org/abs/). It includes both the split conformal inference and cross-validation+. For each type of conformal inference, both conformalized quantile regression (CQR) and standard conformal inference are supported. It provides a pool of convenient learners and allows flexible user-defined learners for conditional mean and quantiles. 

- `conformalCf()` produces intervals for counterfactuals or outcomes with missing values in general.
- `conformalIte()` produces intervals for individual treatment effects with a binary treatment under the potential outcome frameworks. 
- `conformal()` provides a generic framework of weighted conformal inference for continuous outcomes.
- `conformalInt()` provides a generic framework of weighted conformal inference for interval outcomes.

## Installation         

```
# install.packages("devtools")
devtools::install_github("lihualei71/cfcausal")
```
