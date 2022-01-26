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
Save these in a single Rdata file as `Data/Data.Rdata`.

## Running order  

`MarginalAnalysis/` - Scripts in this directory are used to fit the marginal model described in Section 2.1 and 4.2 of the paper <ol>
          <li> `sim_event.R` - Simulate extreme events, see Figure 3 </li>
          <li> `GAM_fit.R` - Marginal GPD, quantile and logistic GAM fits </li>
          <li> `marginal_transform.R` - Transforms data to Laplace margins  </li>
          </ol>

`DependenceAnalysis/` - Scripts in this directory are used to fit the various extremal dependence models described in Section 2.2 and 4.3 <ol>
          <li> `free_fit.R` - Provides the diagnostic "free" fits displayed in Figure 2 </li>
          <li> `spatial_fit.R` - Fit the full spatial model using the censored pseudo-likelihood described in Sections 3.1 and 3.2. Both the full model and the asymptotically dependent model can be fitted </li> 
          <li> `sim_event.R` - Simulate extreme events, see Figure 3 </li>
          </ol>

`AggregateAnalysis/` - Scripts in this directory are used to derive samples of spatial aggregates (denoted $R_\mathcal{A}$ in the paper) and provide fit diagnostics <ol>
          <li> `agg_sim.R` - Draw samples of $R_\mathcal{A}$ </li>
          </ol>

## Uncertainty
To quantify uncertainty in the marginal and extremal dependence fits, replace `Data` with a bootstrap sample of the observations; the bootstrap sample should also be a `n` by `d` matrix. In Richards et al. (2022), we apply the stationary bootstrap (Politis and Romano, 1994) with expected block size corresponding to 48 hours. The function used to derive a single stationary bootstrap sample is given in `src/stat_boot.R`.

