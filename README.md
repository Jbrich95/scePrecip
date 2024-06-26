# Modelling precipitation and spatial aggregate extremes using a conditional framework

Functions for modelling the extremes of spatial aggregates of precipitation using the spatial conditional extremes framework proposed by <a href = "https://www.sciencedirect.com/science/article/pii/S2211675322000471">Wadsworth and Tawn (2022)</a>. 
The provided code is in support of the following two papers:
<ul> 
          <li> Richards, J.,  Tawn, J. A., Brown, S. (2022a). Modelling extremes of spatial aggregates using conditional methods. <i>Annals of Applied Statistics</i>, 16 (4) 2693 - 2713. <u><a href="https://doi.org/10.1214/22-AOAS1609" download>doi.org/10.1214/22-AOAS1609</a></u> </li>
          <li> Richards, J., Tawn, J. A., Brown, S. (2022b). Joint estimation of extreme spatially aggregated precipitation at different scales through mixture modelling. <i>Spatial Statistics</i>, 53:100725. <u><a href="https://doi.org/10.1016/j.spasta.2022.100725" download>doi.org/10.1016/j.spasta.2022.100725</a></u> </li>
</ul>
The latter paper uses a very similar model to the former, but the methodology is applied separately to precipitation classified as being either convective or non-convective; this classification is performed using `conv_identification_algo.R`.


The general framework for fitting the models described in Richards et al. (2022a) and Richards et al. (2022b) is the same for both papers. The differences between the marginal and extremal dependence models are implemented within individual scripts. Following the running order of the scripts, provided below, allows the user to fit the marginal and extremal dependence models to `Data`. Note that if applying the methodology from Richards et al. (2022a), then `Data` is all observations. To apply the mixture model described in Richards et al. (2022b), `Data` will be replaced with either `conv.precip` or `nonconv.precip`, which are outputs from running the script `conv_identification_algo.R`; this is Algorithm 1 in Richards et al. (2022b). 

The algorithm `conv_identification_algo.R` takes in as input - <ul> 
          <li> `Data.grid`: an `M_1` by `M_2` by `n` array of observations. This corresponds to `n` observations on an `M_1` by `M_2` regular grid of spatial locations. </li>
          <li> `lonlat.grid`: an `M_1` by `M_2` by `2` array of lon/lat coordinates. The `[i,j,]`-th element corresponds to the lon/lat coordinates for the location that observes the time series in the `[i,j,]`-th element of `Data.grid`, </li>
</ul>
which is the not the same as the following scripts. The script `conv_identification_algo.R` will output the correct input for the following running order.


## Running order  

Required input - <ul> 
          <li> `Data`: a `n` by `d` matrix of observations. Each row is a time series of length `n` for one of `d` spatial locations. </li>
          <li> `coords`: a `d` by `2` matrix of lon/lat coordinates. The `i`-th row should correspond to the lon/lat coordinates for the location that observes the time series in the `i`-th row of `Data`. </li>
            <li> (if following Richards et al., 2022b) `elev`: a `d` vector of elevation values. The `i`-th element should correspond to the elevation at the `i`-th row of `coords`. </li>
</ul>

Save these in a single Rdata file as `Data/Data.Rdata`. If using the mixture model in Richards et al. (2022b), replace `Data/Data.Rdata` with either `Data/conv.Rdata` or `Data/nonconv.Rdata` and run the code in `MarginalAnalysis/`  and `DependenceAnalysis/` twice - once for each mixture component.

`MarginalAnalysis/` - Scripts in this directory are used to fit the marginal model described in Sections 2.1 and 4.2 of Richards et al. (2022a) or the extensions described in Section 3.2 of Richards et al. (2022b)<ol>
          <li> `GAM_fit.R` - Marginal GPD, quantile and logistic GAM fits </li>
          <li> `marginal_transform.R` - Transforms data to Laplace margins  </li>
          </ol>

`DependenceAnalysis/` - Scripts in this directory are used to fit the various extremal dependence models described in Sections 2.2 and 4.3 of Richards et al. (2022a) or in Section 3.3 of Richards et al. (2022b) <ol>
          <li> `free_fit.R` - Provides the diagnostic "free" fits displayed in Figure 2 in Richards et al. (2022a). Not required if adopting the methods of Richards et al. (2022b)</li>
          <li> `spatial_fit.R` - Fit the full spatial model using the censored pseudo-likelihood described in Sections 3.1 and 3.2 of Richards et al. (2022a). Both the full model and the asymptotically dependent model can be fitted. Also fits convective and nonconvective SCE models from Richards et al. (2022b) </li> 
          <li> `sim_event.R` - Simulate extreme events, see Figure 3 in Richards et al. (2022a) or Figure 6 in Richards et al. (2022b) </li>
          </ol>

`AggregateAnalysis/` - Scripts in this directory are used to derive samples of spatial aggregates (denoted $R_\mathcal{A}$ in the paper) and provide diagnostics and return level estimates <ol>
          <li> `agg_sim.R` - Draw samples of $R_\mathcal{A}$ if using methodology of Richards et al. (2022a); if using the methodology of Richards et al. (2022b), instead run `agg_sim_mix.R` </li> 
           <li> `agg_diags.R` - Produce Q-Q diagnostic plots and estimate return level curves for $R_\mathcal{A}$</li>
          </ol>

## Uncertainty estimates
To quantify uncertainty in the marginal and extremal dependence fits, replace `Data` with a bootstrap sample of the observations; the bootstrap sample should also be a `n` by `d` matrix. In Richards et al. (2022a,2022b), we apply the stationary bootstrap (Politis and Romano, 1994) with expected block size corresponding to 48 hours. The function used to derive a single stationary bootstrap sample is given in `src/stat_boot.R`.

