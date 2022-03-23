# BEST

Bayesian Estimation of the Site frequency specTrum

This software estimates the site-frequency spectrum under a Moran model with biased and reversible mutations and selection using a Bayesian inference. 
 
## The control file 
The control file is where you indicate the location of the data, the model to use, and specify the priors to use during Bayesian inference.


| Parameter | Description |
|---|---|
| ```model``` | ```S``` stands for the Moran model with reversible mutations and selection whereas ```N``` stands for the neutral Moran model. |
| ```site_frequency_spectrum_file``` | The location of the sampled site frequency spectrum. |
|```mutation_rate_mean```| The mutation rate. We advise to use mutation rate estimated independently of previous estimate of the effective population size (e.g., pedigree-based, etc.). The problem of using such mutation rates is that the inferences will be circular, as we are jointly estimating the effective population size in our analyses. |
|```mutation_rate_variance``` | The variance of the mutation rate. It should be given from the place you have taken the mutation rate. But if not, you can came up with a suitable range. In order to facilitate the calculations, the prior on the mutation range is normal distributed. So you can set this parameters according to your expectations.|
|```population_size_dgamma_prior_mean```| The mean of the discrete prior gamma distribtuion set to the population size. A good option for a uninformative gamma priors is setting alpha=1000 and beta=1000. This is equivalent to set the mean of the gamma distribtuion to 1. The mean and the variance of a discrete gamma distribtuion relate to alpha and beta via $\alpha = E*E/V$ and $\beta = E/V$ (where E stands for the mean and V the variance).  |
|```population_size_dgamma_prior_variance```| The variance of the discrete prior gamma distribtuion set to the population size. A typical good option for uninformative priors with gamma distribution is setting alpha=1000 and beta=1000. This is equivalent to set the variance of the gamma distribtuion to 1000. |
|```fitness_coefficient_gamma_prior_mean```| The mean of the prior gamma distribution set to the fitness coefficient. |
|```fitness_coefficient_gamma_prior_variance``` | The variance of the prior gamma distribtuion set to the fitness coefficient. |
|```number_generations```| The number of generations for which the analysis will be run.|
|```number_chains```| The number of independent runs  starting from different parameter set. |
|```sample_frequency```| The sample frequency determines how often the chain is sampled.|


## The sampled site frequency spectrum 
This software uses the frequency spectrum of two variants to produce inferences. Let us imagine that we have a determining number of sequences aligned for a given number of individuals *S*. For each genomic position, the number of times the variant *A* is present among the sampled individuals can be counted. The distribution of these frequencies is the site frequency spectrum. This frequency can vary between 0 (when none of the individuals have the variant *A*) and *S* (when all the individuals have the *A* variant). A similar reasoning can be made for variant *B*.

If the individuals are diploid then we count chromosomes, so the sample size will be *2S*; and the same follows for other ploidies. These software only needs a vector of counts separated by spaces. Let us consider the following example:
```
10 8 5 6 9
````
From this count vector, the  immediately assumes that a total of 38 genomic positions were observed: 10 of which where the variant *A* was not observed among the individuals (or chromosomes), and 9 where all individuals had variant *A*. As the vector has 5 elements, it assumes that four chromosomes were counted, i.e., there are 4 haploid individuals in the sample. 

## Compiling and Running BEST

Now that the control file and the sampled site frequency spectrum are created, let us compile and runs BEST.


## The output file

## Cite us!

In preparation
