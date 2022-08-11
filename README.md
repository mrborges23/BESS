# BESS

Bayesian Estimator of the Site-frequency Spectrum

A bayesian estimator of the site-frequency spectrum under a Moran model with biased and reversible mutations and selection. 

 

## Version 

0.1 (July 2022)



## Citation

*Title* Borges, Kotari, Bergman, Chase, Mugal and Kosiol Birxiv DOI XXX




## Downloading and compiling BESS

First, we need to download all the necessary files to compile and run BESS, which are in this GitHub repository. You can download these files manually or directly from the terminal:

```

git clone https://github.com/mrborges23/bess.git bess

```

BESS is implemented in C++ and needs to be compiled. This step will produce an executable that can then be used to run BESS. We use the `g++` compiler, but others could be used as well:

```

g++ bess.cpp -o bess -O3 -std=c++20

```

This step should have created an executable called `bess` in your working directory.

## Control file

The control file includes all the necessary parameters to run BESS: it is where you indicate the location of the data, the model to use, and specify the priors for Bayesian estimation.

* ```model```: ```S``` stands for the Moran model with reversible mutations and selection, whereas ```N``` stands for the neutral Moran model.
* ```site_frequency_spectrum_file```: The location of the sampled site frequency spectrum.
* ```mutation_rate```: The mutation rate. We advise using a mutation rate estimated independently of any previous estimates of effective population size (e.g., pedigree-based, etc.). The problem with using effective-population-size estimated mutation rates, even if obtained from other software, is that BESS jointly estimates the effective population size, making the inferences circular.
* ```population_size_prior_mean```: The prior mean of a normal distribution. We note that each new proposed effective population size during the MCMC is corrected so it cannot take negative or decimal numbers. 
* ```population_size_prior_variance```: The prior variance of a normal distribtuion.
* ```fitness_coefficient_prior_mean```: The prior mean of the prior normal distribution. As for the effective population size, each newly proposed fitness value during the MCMC cannot take negative values.
* ```fitness_coefficient_prior_variance```: The prior variance of the prior normal distribution.
* ```number_generations```: The number of generations for which the MCMC is intended to run.
* ```number_chains```: The number of independent runs, each starting from a different parameter set.
* ```sample_frequency```: The sample frequency determines how often the chain is sampled.
* ```continuous_approximation```: This option tells BESS if a continuous approximation of the site-frequency spectrum should be employed and how many points should be used in such approximation. This approximation can be useful with the model with selection with populations of larger effective population size, for which the likelihood computation takes considerable time. If this value differs from 0, the approximation is employed with the specified number of points. We recommend using 10 000 points, which is already a good fit for the site-frequency spectrum.


## The sampled site frequency spectrum 

This software uses the site frequency spectrum to infer several population parameters. But what is the site frequency spectrum? Let us imagine that we have sequenced *S* genomic positions from a population of *I* individuals. Each genomic position can only have the allele *A* or *B*. If we count the number of times the variant *A* is present among the sampled individuals for each genomic position. Such a count can vary between 0 (none of the individuals have the variant *A*) and *S* (all the individuals have the variant *A*). The frequency (among the genomic positions) of the different allelic counts (among the individuals) constitutes the (sampled) site frequency spectrum. An example will clarify this hideous sentence. Let us consider the following toy site-frequency spectrum:

```

10 8 5 6 9

````
We know immediately that a total of 38 genomic positions were sequenced (by summing all the counts): 10 of which the variant *A* was not observed, 9 of which all individuals had the variant *A*. As the vector has five elements, it tells us that four individuals were counted, i.e., there are four haploid individuals in the sample. We note that if the individuals are diploid, then we count chromosomes instead of individuals, so the sample size will be *2S*. The same rationale follows for other ploidies. 

BESS requires a vector of counts separated by spaces, which location is given in the control file (```site_frequency_spectrum_file```). 



## Running BESS

To run BESS, simply place your sampled site-frequency spectrum and the control file together with the BESS executable. Then, open the terminal and run the executable `bess` followed by the name of the control file:

```

./bess control.cf

```

BESS immediately prints out a short description of the data file. Please, confirm that this information conforms to your data. To ensure BESS is running, you should get the message `BESS has started!`. When the analyses finish, you should get the message `BESS has finished!`. BESS periodically writes to the `output_file`.

## Output file

The output file (```.log``` file) includes information on the estimated parameters. Each line is a sample from the posterior. 

```
gen	  lnL	         muji	       muij	         sigmaj	       N
100	  -7.52836e+06	 2.3684e-08	   2.64722e-08	 6.77759e-08   13710
200	  -7.48412e+06	 2.38531e-08   2.62628e-08	 6.77759e-08   1015
300	  -7.41649e+06	 2.41379e-08   2.5926e-08	 6.77759e-08   1015
400	  -7.33793e+06	 2.44865e-08   2.55355e-08	 6.77759e-08   1015
500	  -7.20977e+06	 2.51017e-08   2.48991e-08	 8.61609e-09   12724
```

To estimate the population parameters, we need first to check the convergence and mixing of MCMC chains. In this GitHub repository, you find a site-frequency spectrum and a control file that should produce output files which trace plots should be similar to the following ones:

![](https://github.com/mrborges23/BESS/blob/main/mcmc_plots.png)

These are examples of well converged and mixed chains. To estimate each parameter, the burn-in phase, in which the parameters quickly change their values, should be first removed. Removing the first 10% generations would be enough to eliminate the burn-in phase for these runs, but thei t should be checked for each data set, as it may vary considerably. 

Note: This data comes from a population of 11 Sumatra-Orang-Utans; as they are diploid, 23 the site-frequency has 23 entries. 

## Compilation errors and how to solve them

So far, we have found none but let us know otherwise.


## Questions and bug reporting

Please use **Issues** to report possible bugs, suggest enhancement features, or if you need help using BESS. If you have more theoretical or biological questions, you can directly contact Rui Borges (ruiborges23@gmail.com).

We noticed that some runs get stuck during the MCMC; this is noticeable because the parameters do not change. These are related to numerical issues that tend to happen when the site-frequency spectrum has few individuals (two or three individuals). We regret not having a fancier solution to this problem, but in the situations we encountered, we could always work around this numerical issue by just re-running BESS. 


## License

This program is free software. You can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software. See the GNU General Public License (http://www.gnu.org/licenses/) for more details.
