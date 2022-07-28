# BESS

Bayesian Estimator of the Site-frequency Spectrum

Bayesian estimation of the site-frequency spectrum under a Moran model with biased and reversible mutations, and selection. 

 

## Version 

0.1 (July 2022)



## Citation

In preparation





## Downloading and compiling BESS

First of all, we need to download all the necessary files to compile and run BESS; these files are in this GitHub repository. You can download these files manually or directly from the terminal:

```

git clone https://github.com/mrborges23/best.git best

```

BESS is implemented in C++ and needs to be compiled; this step will produce an executable that can then be used to run BESS. We use the `g++` compiler, but others could be used instead:

```

g++ best.cpp -o best -O3 -std=c++20

```

This step should have created an executable called `best` in your working directory.

## Control file

The control file includes all the necessary parameters to run BESS:

The control file is where you indicate the location of the data, the model to use, and specify the priors to use during Bayesian inference.

* ```model```: ```S``` stands for the Moran model with reversible mutations and selection whereas ```N``` stands for the neutral Moran model.
* ```site_frequency_spectrum_file```: The location of the sampled site frequency spectrum.
* ```mutation_rate```: The mutation rate. We advise using a mutation rate estimated independently of any previous estimate of the effective population size (e.g., pedigree-based, etc.). The problem with using such mutation rates is that the inferences will be circular. Note that Best jointly estimates the effective population size.
* ```population_size_prior_mean```: The prior mean of a normal distribution A good option for a uninformative gamma priors is setting the mean to ... and the variance to ... . The effective population size can only take integer values; we renormalize the normal distribution so that only positive integers are possible to be observed.
* ```population_size_prior_variance```: The prior variance of a normal distribtuion.
* ```fitness_coefficient_prior_mean```: The prior mean of the prior normal distribution. The fiftness coefficients can only take positive values; we renomalized the normal distribution so that it  are possible to be observed.
* ```fitness_coefficient_prior_variance```: The prior variance of the prior normal distribtuion.
* ```number_generations```: The number of generations for which the analysis will be run.
* ```number_chains```: The number of independent runs  starting from different parameter set.
* ```sample_frequency```: The sample frequency determines how often the chain is sampled.
* ```lower_bound```: The sample frequency determines how often the chain is sampled.


## The sampled site frequency spectrum 

This software uses the site frequency spectrum to infer the population parameters. But what is the site frequency spectrum? Let us imagine that we have *S* sequences for different individuals.  For each genomic position, the number of times the variant *A* is present among the sampled individuals can be counted. The distribution of these frequencies among the genomic position is the (sampled) site frequency spectrum. This frequency can vary between 0 (when none of the individuals have the variant *A*) and *S* (when all the individuals have the variant *A*).

If the individuals are diploid, then we count chromosomes, so the sample size will be instead *2S*; and the same follows for other ploidies. BESS software requires a vector of counts separated by spaces, which location is given in the control file (```site_frequency_spectrum_file```). Let us consider the following example:

```

10 8 5 6 9

````

From this count vector, one immediately knows that a total of 38 genomic positions were sequenced: in 10 of which the variant *A* was not observed, and in 9 all individuals had variant *A*. As the vector has 5 elements, it assumes that four individuals (or chomosomes for diploids and other ploidies) were counted, i.e., there are 4 haploid individuals in the sample. 

## Running BESS

To run BESS, simply place your sampled site frequency spectrum and the control file together with the BESS executable. Then, open the terminal and run the executable `best` followed by the name of the control file:

```

./best control.cf

```

BESS immediately prints out a short description of the data file. Please, confirm that this information conforms to your data. To make sure BESS is running, you should get the message `BESS has started!`. When the analyses finish, you should get the message `BESS has finished!`. BESS periodically writes to the `output_file`.

## Output file

The output file (```.log``` file) includes information on the estimated parameters. Each line is a sample from the posterior. Please make sure to exclude some of the initial samples from your estimates. These correspond to the burni phase, where the chain  

```
gen    lnL             muji           muij           N
100    -3.988e+07      8.97193e-09    5.7387e-09     9089
200    -3.98104e+07    8.97649e-09    5.73683e-09    22303
300    -3.97732e+07    8.97649e-09    5.73683e-09    40726
```

An example:

## Compilation errors and how to solve them

So far, we have found none. But let us know if you are experiencing some issues and if so how did you solve them. 


## Questions and bug reporting

Please use **Issues** to report possible bugs, suggest enhancement features, or if you need help using BESS. If you have more theoretical or biological questions, you can directly contact Rui Borges (ruiborges23@gmail.com).

We noticed that some runs get stucked during the MCMC (this is noticeable because the parameters do not change at all; your trace plots would look like this:


These are related to numerical issues that tend to hapen with site frequency spectrum from two or three individuals. We were able to solve these by re-runing BESS. 


## License

This program is free software. You can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software. See the GNU General Public License (http://www.gnu.org/licenses/) for more details.
