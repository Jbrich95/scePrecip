# Modelling of extreme precipitation using spatial conditional extremes
Functions for modelling the extremes of spatial aggregates of precipitation using the spatial conditional extremes framework. The provided code is in support of the paper:
Richards, J., Tawn, J. A., Brown, S. (2022). <i>Modelling Extremes of Spatial Aggregates of Precipitation using Conditional Methods</i>, <a href = "https://arxiv.org/pdf/2102.10906.pdf">ArXiv.</a>
## Installation

```r
library(devtools)
install_github("https://github.com/Jbrich95/scePrecip")
```

## Running order
Required input - `Data`: a $`n \times d`$ matrix of observations. 
<ol>
  <li>MarginalAnalysis/ - Scripts in this directory are used to fit the marginal model described in Section 2.1 and 4.2 of the paper.<ol>
<li>GAM_fit.R - Marginal GPD, quantile and logistic GAM fits</li>
<li>marginal_transform.R - Transforms data to Laplace margins</li></ol>
</li>
  <li>DependenceAnalysis/ - Scripts in this directory are used to fit the various extremal dependence models described in Section 2.2 and 4.3 of the paper<ol>
<li>free_fit.R - Provides the diagnostic "free" fits displayed in Figure 2</li>
<li>spatial_fit_AI.R - Fit the full spatial asymptotic model using the censored pseudo-likelihood described in Sections 3.1 and 3.2 </li></ol>
</li>
  
</ol>
