# wispack
 Rcpp implementation of warped-sigmoid Poisson-process mixed-effect models

# idea of pinned-age test: to test whether model is underconstrained, pin ages to the values found for them 
  in the full fit while bootstrapping. Presumably the optimal fits on the resamples will have
  different values. If the model is underconstrained and can make up for these sub-optimal values 
  by adjusting other parameters (presumably the random effects), then the bootstrapped fits will be
  just as good, whether or not the ages are pinned. 
  
# How it works: 
  - age_effect_mask covers the fitted parameter vector.
  - the boundary penalty function holds any parameters exposed as "true" by this mask to a value 
    set in pinned_ages.
  - for the initial fit, the age_effect_mask is all false, so there is no effect.
  - then, all age-related parameters are saved from that initial optimization into the pinned_ages 
    vector and set to "true" in the age_effect_mask.
  - So, as the bootstraps are fit, all age-related parameters are held to the values found in the initial fit, 
    via penalization for moving away from those values. 