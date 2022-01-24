# Modelling of extreme precipitation using spatial conditional extremes
Functions for modelling the extremes of spatial aggregates of precipitation using the spatial conditional extremes framework. The provided code is in support of the paper:
Richards, J., Tawn, J. A., Brown, S.(2022)  <i>Modelling Extremes of Spatial Aggregates of Precipitation using Conditional Methods</i>, <a href = "https://arxiv.org/pdf/2102.10906.pdf">ArXiv.</a>
## Installation

```r
library(devtools)
install_github("https://github.com/Jbrich95/scePrecip")
```

## Running order
<ol>
  <li>MarginalAnalysis/<ol>
<li>GAM_Fit.R</li>
<li>Marg_Transform.R</li></ol>
</li>
  <li>DependenceAnalysis/<ol>
<li>Free_fits.R</li>
<li>Marg_Transform.R</li></ol>
</li>
  
</ol>
