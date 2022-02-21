# SplitKnockoff
This is a R package for Split Knockoff method, related paper can be found as below.
>Yang Cao, Xinwei Sun and Yuan Yao, Controlling the False Discovery Rate in Structural Sparsity: Split Knockoffs, [arXiv:2103.16159](https://arxiv.org/abs/2103.16159)

# Vignette for this package
Please visit website for vivid illustration [Vignette](http://www.wanghaoxue.com/4/splitknockoff.html)

## For Installation 
```R
install.packages("SplitKnockoff")
library(SplitKnockoff)
```
Reference manual for this package can be found [SplitKnockoff.pdf](https://github.com/wanghaoxue0/SplitKnockoff/blob/main/SplitKnockoff_0.1.0.pdf)

Split Knockoff is a novel method for controlling the false discovery rate (FDR) in structural sparsity setting. This proposed scheme relaxes the linear subspace constraint to its neighborhood, often known as variable splitting in optimization. **This vignette illustrates the usage of Splitknockoff with simulation experiments and will help you apply the split knockoff method in a light way.**

```R
install.packages("SplitKnockoff")   # just one line code to install our package
```

There is another Matlab package available. Simulation experiment of both packages verify the conclusion of the paper. This R package is more convenient when applied to new dataset as **glmnet** can be used directly by 

```R
install.packages("glmnet'")  # just one line code to install the glmnet tool
```

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




