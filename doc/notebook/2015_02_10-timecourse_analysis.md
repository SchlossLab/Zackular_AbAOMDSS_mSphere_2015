# Time course analysis


**Purpose:** We want to characterize the change in the microbiome over the
course of the model


### Setup


```r
library(randomForest, quietly=TRUE)
library(knitr, quietly=TRUE)
sessionInfo()

set.seed(6201976)

source("code/rf_baseline_analysis.R")
```

```
## R version 3.1.2 (2014-10-31)
## Platform: x86_64-apple-darwin13.4.0 (64-bit)
## 
## locale:
## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
## 
## attached base packages:
## [1] stats     graphics  grDevices utils     datasets  methods   base     
## 
## other attached packages:
## [1] randomForest_4.6-10 knitr_1.8          
## 
## loaded via a namespace (and not attached):
## [1] evaluate_0.5.5 formatR_1.0    stringr_0.6.2  tools_3.1.2
```



First, we want to read in the distance matrix and see which treatments are
changing the most...


```
##      Group.1   x.Min. x.1st Qu. x.Median   x.Mean x.3rd Qu.   x.Max.
## 1     AllAbs 0.001237  0.007671 0.030260 0.208200  0.458500 0.614700
## 2      Metro 0.489800  0.588500 0.641800 0.619900  0.646700 0.732600
## 3      NoAbs 0.326700  0.571600 0.669600 0.654000  0.780500 0.914500
## 4      Strep 0.110600  0.166600 0.174000 0.331900  0.538700 0.669700
## 5 StrepMetro 0.679100  0.742100 0.795900 0.779900  0.827400 0.855100
## 6       Vanc 0.070500  0.082360 0.094210 0.116000  0.138700 0.183200
## 7  VancMetro 0.068550  0.103500 0.138500 0.130400  0.161400 0.184300
## 8  VancStrep 0.857500  0.872500 0.888200 0.903900  0.936600 0.964600
```

Was surprised that the "AllAbs" seemed bimodal. Of the mice, 5 didn't move at
all and 3 moved a fair amount. Curious to see how well they replicated each other
at the beginning and the end of the model.



```
##           AllAbs  Metro  NoAbs   Strep StrepMetro    Vanc VancMetro
## Min.    0.001354 0.1307 0.1569 0.07413   0.008213 0.03332   0.03032
## 1st Qu. 0.008829 0.1959 0.2598 0.09710   0.058240 0.03365   0.03979
## Median  0.015120 0.2389 0.4541 0.12920   0.152600 0.03398   0.04926
## Mean    0.055430 0.2332 0.4285 0.15360   0.187100 0.04207   0.06302
## 3rd Qu. 0.111500 0.2520 0.5480 0.21040   0.300600 0.04644   0.07938
## Max.    0.161300 0.3549 0.7683 0.28080   0.472900 0.05890   0.10950
##         VancStrep
## Min.      0.07723
## 1st Qu.   0.18610
## Median    0.38480
## Mean      0.44140
## 3rd Qu.   0.71920
## Max.      0.81380
```

Again, I was struck by the variation among the Day 0 samples for all
antibiotics. Let's see what it looks like at the end of the model...


```
##           AllAbs  Metro  NoAbs   Strep StrepMetro   Vanc VancMetro
## Min.    0.001725 0.1143 0.2334 0.03205    0.01524 0.1006   0.01764
## 1st Qu. 0.031450 0.1822 0.5819 0.16300    0.04075 0.1067   0.07927
## Median  0.448100 0.2883 0.6930 0.30720    0.06168 0.1128   0.17660
## Mean    0.324300 0.2899 0.6695 0.30440    0.12850 0.1580   0.20160
## 3rd Qu. 0.596200 0.3315 0.7928 0.42700    0.24280 0.1867   0.33230
## Max.    0.665000 0.6047 0.9250 0.58280    0.28290 0.2605   0.40790
##         VancStrep
## Min.      0.02841
## 1st Qu.   0.04940
## Median    0.06062
## Mean      0.07115
## 3rd Qu.   0.07190
## Max.      0.16450
```

More variation at the end point.

What happens when we plot the distance changed between day 0 and the end with
the number of tumors that were observed?



```
## Error in `[.data.frame`(tumor_counts, names(start_end_dist)): undefined columns selected
```

<img src="results/figures/distance_tumor_corr-1.png" title="plot of chunk distance_tumor_corr" alt="plot of chunk distance_tumor_corr" style="display: block; margin: auto;" />

My eye wants to see some type of saturation curve where the more distance
between the start and end point, the more tumors occur. But... the $\rho$
statistic was 0.26 (P=
0.08).

Let's break down and build an ordination...

<img src="results/figures/unnamed-chunk-1-1.png" title="plot of chunk unnamed-chunk-1" alt="plot of chunk unnamed-chunk-1" style="display: block; margin: auto;" />

Let's see whether we can build a random forest to predict whether a sample comes
from the beginning or the end of the model.



Looks pretty crappy. Here's the confusion matrix:



|      | begin| end| class.error|
|:-----|-----:|---:|-----------:|
|begin |    35|  11|   0.2391304|
|end   |     8|  40|   0.1666667|



The out of bag error rate was 20.2%.


How about if we do it based on the final samples and see whether there was
overlap with the OTUs we detected at the beginning of the model.

<img src="results/figures/unnamed-chunk-3-1.png" title="plot of chunk unnamed-chunk-3" alt="plot of chunk unnamed-chunk-3" style="display: block; margin: auto;" /><img src="results/figures/unnamed-chunk-3-2.png" title="plot of chunk unnamed-chunk-3" alt="plot of chunk unnamed-chunk-3" style="display: block; margin: auto;" /><img src="results/figures/unnamed-chunk-3-3.png" title="plot of chunk unnamed-chunk-3" alt="plot of chunk unnamed-chunk-3" style="display: block; margin: auto;" />

This model does a pretty good job of predicting the number of tumors, but it
is pretty redundant with the mBio paper and the baseline model.


Let's try building a model to connect change in relative abundance between
baseline and the end point with the number of tumors




This results in 4 OTUs - Ecoli, Bacteroidetes, Odoribacter, and
Beta-proteobacteria. Not sure what it really says...
