# README for the toy model

See draft article for description of the model and methods

The main script is 'naiveEstimationExp.m'. 

Running 'naiveEstimationExp.m' will:
* Define the true mean of a and b
* Draw individual values for a and b for a specified number of individuals from the population level distribution, assumed to be independent exponentials
* Generate synthetic data by solving the forward model for each individual and adding IID Gaussian noise on the log scale
  
* Define a prior for a and b (currently uniform betwene specified bounds)
* Generate a specified number of samples from the prior
* For each sample, evaluate the Monte Carlo likelihood function (see Methods)
* Plot a colour map on the likelihood values in bivariate parameter space ( mean(a) , mean(b) )



Some things on the to do list:
* Investigate the shape of likelihood surface for diffferent numbers of individuals (10, 100, 1000?)
* Investigate the shape of likelihood surface for different values of the ratio mean(a)/mean(b), with the hypothesis being that ratios close to will be non-identifiable but ratios away from 1 may be more identifiable.

* Run the model in monolix? 


