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


```r
opts_chunk$set(dev = c("png", "cairo_pdf"))
opts_chunk$set(results = "hold")
opts_chunk$set(fig.show = "hold")
opts_chunk$set(warning = FALSE)
opts_chunk$set(fig.align = "center")
opts_chunk$set(echo = FALSE)
opts_chunk$set(cache = FALSE)
```
