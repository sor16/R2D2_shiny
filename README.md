## R2-D2 prior distribution visualization
The R2-D2 prior, presented by Zhang et al, provides an elegant way to define a prior on regression coefficients jointly by supplying a prior on $R^2$, the coefficient of variation, and then distributing through to the coefficient through a Dirichlet distribution. In addition to being an intuitive way to define a prior in a regression setting, an appealing property of the R2-D2 is that it can provide shrinkage. The aim of this Shiny app is to visualize how different choices of parameters to the prior distribution affects the resulting prior on the regression coefficients in an interactive way.

## References
Yan Dora Zhang, Brian P. Naughton, Howard D. Bondell & Brian J. Reich (2022) Bayesian Regression Using a Prior on the Model Fit: The R2-D2 Shrinkage Prior, Journal of the American Statistical Association, 117:538, 862-874, DOI: 10.1080/01621459.2020.1825449
