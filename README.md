# SplitKnockoff
#### Introduction

Split Knockoff is a data adaptive variable selection framework for controlling the (directional) false discovery rate (FDR) in structural sparsity, where variable selection on linear transformation of parameters is of concern. This proposed scheme relaxes the linear subspace constraint to its neighborhood, often known as variable splitting in optimization, which leads to better incoherence and possible performance improvement in both power and FDR. **This vignette illustrates the usage of Splitknockoff with simulation experiments and will help you apply the split knockoff method in a light way.** On top of that, **an application to Alzheimer’s Disease** study with MRI data is given as an example to illustrate that the split knockoff method can disclose important lesion regions in brains associated with the disease and connections between neighboring regions of high contrast variations during disease progression.

```R
install.packages("SplitKnockoff")   # just one line code to install our package
```

This is a R implement on the Matlab version of Split Knockoffs. This R package is more convenient as **glmnet** can be used directly by 

```R
install.packages("glmnet'")  # just one line code to install the glmnet tool
```

**Please update your glmnet package to the latest version for smooth usage of this package**, the examples illustrated here are tested with *glmnet 4.2*.

For more information, please see the manual inside this package.

#### Key function

**sk.filter(X, D, y, option)**   : the main function, Split Knockoff filter, for variable selection in structural sparsity problem.

#### function involved frequently

**sk.create(X, y, D, nu, option)**: generate the split knockoff copy for structural sparsity problem.

**W_sign(X, D, y, nu, option)**: generate the W statistics for split knockoff on a split LASSO path.

**sk.select(W, q)**: this function is for variable selection based on W statistics.



**Please update your glmnet package to the latest version for smooth usage of this package**, this simulation uses *glmnet 4.2*





