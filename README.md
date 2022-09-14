# Modelling extreme precipitation and spatial aggregates using spatial conditional extremes
Functions for modelling the extremes of spatial aggregates of precipitation using the spatial conditional extremes framework. The provided code is in support of the papers:
<ul> 
          <li> Richards, J., Tawn, J. A., Brown, S. (2022a). <i>Modelling Extremes of Spatial Aggregates of Precipitation using Conditional Methods</i>, <a href = "https://arxiv.org/pdf/2102.10906.pdf">arXiv.</a> </li>
          <li> Richards, J., Tawn, J. A., Brown, S. (2022b). <i>Joint Estimation of Extreme Spatially Aggregated Precipitation at Different Scales through Mixture Modelling</i>, <a href = "https://arxiv.org/pdf/2111.08469.pdf">arXiv.</a> </li>
</ul>
The latter paper uses a very similar model, but it is applied separately to precipitation classified as being either convective or non-convective; this is performed using `conv_identification_algo.R`.

The general framework for fitting the models described in Richards et al. (2022a) and Richards et al. (2022b) is the same for both papers. The differences between the marginal and extremal dependence models are implemented within individual scripts. Following the running order of the scripts, provided below, allows the user to fit the marginal and extremal dependence models to `Data`. Note that if applying the methodology from Richards et al. (2022a), then `Data` is all observations. To apply the mixture model described in Richards et al. (2022b), replace `Data` with either `conv.precip` or `nonconv.precip`, which are outputs from running the script `conv_identification_algo.R`; this is Algorithm 1 in Richards et al. (2022b). 

The algorithm `conv_identification_algo.R` takes in as input:
<ul> 
          <li> `Data.mat`: an `n` by `M_1` by `M_2` array of observations. This corresponds to `n` observations on an `M_1` by `M_2` regular grid of spatial locations. </li>
          <li> `lonlat.mat`: an `M_1` by `M_2` by `2` array of lon/lat coordinates. The `[i,j,]`-th element corresponds to the lon/lat coordinates for the location that observes the time series in the `[,i,j]`-th element of `Data.mat`, </li>
</ul>
which is the not the same as the following scripts. The script `conv_identification_algo.R` will output the correct input for the following running order.


## Running order  

Required input - <ul> 
          <li> `Data`: a `n` by `d` matrix of observations. Each row is a time series of length `n` for one of `d` spatial locations. </li>
          <li> `coords`: a `d` by `2` matrix of lon/lat coordinates. The `i`-th row should correspond to the lon/lat coordinates for the location that observes the time series in the `i`-th row of `Data`. </li>
</ul>

Save these in a single Rdata file as `Data/Data.Rdata`. If using the mixture model in Richards et al. (2022b), replace `Data/Data.Rdata` with either `Data/conv.Rdata` or `Data/nonconv.Rdata`.

`MarginalAnalysis/` - Scripts in this directory are used to fit the marginal model described in Sections 2.1 and 4.2 of Richards et al. (2022a) or the extensions described in Section 3.2 of Richards et al. (2022b)<ol>
          <li> `GAM_fit.R` - Marginal GPD, quantile and logistic GAM fits </li>
          <li> `marginal_transform.R` - Transforms data to Laplace margins  </li>
          </ol>

`DependenceAnalysis/` - Scripts in this directory are used to fit the various extremal dependence models described in Sections 2.2 and 4.3 of Richards et al. (2022a) or in Section 3.3 of Richards et al. (2022b <ol>
          <li> `free_fit.R` - Provides the diagnostic "free" fits displayed in Figure 2 in Richards et al. (2022a)</li>
          <li> `spatial_fit.R` - Fit the full spatial model using the censored pseudo-likelihood described in Sections 3.1 and 3.2 of Richards et al. (2022a). Both the full model and the asymptotically dependent model can be fitted </li> 
          <li> `sim_event.R` - Simulate extreme events, see Figure 3 in Richards et al. (2022a)</li>
          </ol>

`AggregateAnalysis/` - Scripts in this directory are used to derive samples of spatial aggregates (denoted $R_\mathcal{A}$ in the paper) and provide diagnostics and return level estimates <ol>
          <li> `agg_sim.R` - Draw samples of $R_\mathcal{A}$ if using methodology of Richards et al. (2022a); if using the methodology of Richards et al. (2022b), this will instead give the corresponding contribution to $R_\mathcal{A}$ of either `conv.precip` or `nonconv.precip` (see Figure S9 in the supplement), depending on the choice of `Data` </li> 
          <li> If following Richards et al. (2022b), re-run all scripts starting from `GAM_fit.R` using whichever of `Data=conv.precip` or `Data=nonconv.precip` that was not used initially.</li>
           <li> `agg_diags.R` - Produce Q-Q diagnostic plots and estimate return level curves for $R_\mathcal{A}$</li>
          </ol>

## Uncertainty estimates
To quantify uncertainty in the marginal and extremal dependence fits, replace `Data` with a bootstrap sample of the observations; the bootstrap sample should also be a `n` by `d` matrix. In Richards et al. (2022a,2022b), we apply the stationary bootstrap (Politis and Romano, 1994) with expected block size corresponding to 48 hours. The function used to derive a single stationary bootstrap sample is given in `src/stat_boot.R`.

