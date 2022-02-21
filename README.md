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

#### Simulation Results

Figure 1b of the paper, comparing the different distributions of Knockoff statistics and Split Knockoff statistics.

<img src="https://tva1.sinaimg.cn/large/008i3skNgy1gtcsdcbx6aj60rs0rs0tp02.jpg" alt="figure_1b" style="zoom:30%;" />

Figure 2a of the paper, **FDR** in **$D_1$**

<img src="https://tva1.sinaimg.cn/large/008i3skNgy1gtd5rf6auuj60rs0rsdgw02.jpg" alt="figure_11" style="zoom:30%;" />

Figure 2b of the paper, **Power** in **$D_1$**

<img src="https://tva1.sinaimg.cn/large/008i3skNgy1gtd5rm70nnj60rs0rs75g02.jpg" alt="figure_12" style="zoom:30%;" />

Figure 2b of the paper, **FDR** in **$D_2$**

<img src="https://tva1.sinaimg.cn/large/008i3skNgy1gtd5ru8jkyj60rs0rsab702.jpg" alt="figure_21" style="zoom:30%;" />

Figure 2b of the paper, **Power** in **$D_2$**

<img src="https://tva1.sinaimg.cn/large/008i3skNgy1gtd5s6tpfrj60rs0rsmyi02.jpg" alt="figure_22" style="zoom:30%;" />

Figure 2b of the paper, **FDR** in **$D_3$**

<img src="https://tva1.sinaimg.cn/large/008i3skNgy1gtd5sfaxpyj60rs0rs0tr02.jpg" alt="figure_31" style="zoom:30%;" />

Figure 2b of the paper, **Power** in **$D_3$**

<img src="https://tva1.sinaimg.cn/large/008i3skNgy1gtd5sma46bj60rs0rs0u502.jpg" alt="figure_32" style="zoom:30%;" />

Figure 3a of the paper,  *Effect of Feature Correlation:  FDR for* q = 0.2

<img src="https://tva1.sinaimg.cn/large/008i3skNgy1gtcshfwrffj60rs0rsq4402.jpg" alt="figure_3a" style="zoom:30%;" />

Figure 3b of the paper,  *Effect of Feature Correlation: Power for* q = 0.2

<img src="https://tva1.sinaimg.cn/large/008i3skNgy1gtcshwgh4sj60rs0rsgn202.jpg" alt="figure_3b" style="zoom:30%;" />

Figure 4a of the paper,  *Effect of SNR:  FDR for* q = 0.2

<img src="https://tva1.sinaimg.cn/large/008i3skNgy1gtcskl7xfrj60rs0rsmy602.jpg" alt="figure_4a" style="zoom:30%;" />

Figure 4b of the paper,  *Effect of SNR: Power for* q = 0.2

<img src="https://tva1.sinaimg.cn/large/008i3skNgy1gtcsmudg78j60rs0rsdgy02.jpg" alt="figure_4b" style="zoom:30%;" />

Figure 5a of the paper, *Effect of Sparsity Level: FDR for* q = 0.2

<img src="https://tva1.sinaimg.cn/large/008i3skNgy1gtcso1dmm2j60rs0rs3zm02.jpg" alt="figure_5a" style="zoom:30%;" />

Figure 5b of the paper, *Effect of Sparsity Level: Power and FDR for* q = 0.2

<img src="https://tva1.sinaimg.cn/large/008i3skNgy1gtcspcan83j60rs0rsabf02.jpg" alt="figure_5b" style="zoom:30%;" />




