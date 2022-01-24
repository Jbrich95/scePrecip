# Modelling of extreme precipitation using spatial conditional extremes
Functions for modelling the extremes of spatial aggregates of precipitation using the spatial conditional extremes framework. The provided code is in support of the paper:
Richards, J., Tawn, J. A., Brown, S. (2022). <i>Modelling Extremes of Spatial Aggregates of Precipitation using Conditional Methods</i>, <a href = "https://arxiv.org/pdf/2102.10906.pdf">ArXiv.</a>
## Installation

```r
library(devtools)
install_github("https://github.com/Jbrich95/scePrecip")
```
Required input - <ul> 
          <li> `Data`: a `n` by `d` matrix of observations. Each row is a time series of length `n` for one of `d` spatial locations. </li>
          <li> `coords`: a `d` by `2` matrix of lon/lat coordinates. The `i`-th row should correspond to the lon/lat coordinates for the location that observes the time series in the `i`-th row of `Data`. </li>
</ul>

## Running order     

<ol>
         <li> MarginalAnalysis/ - Scripts in this directory are used to fit the marginal model described in Section 2.1 and 4.2 of the paper.
                    <ol>
                              <li> GAM_fit.R - Marginal GPD, quantile and logistic GAM fits </li>
                              <li> marginal_transform.R - Transforms data to Laplace margins </li>
                    </ol>
          </li>
          
      
    
</ol>

end
