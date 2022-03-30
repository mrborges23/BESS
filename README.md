# BEST

Bayesian Estimation of the Site frequency specTrum

This software estimates the site-frequency spectrum under a Moran model with biased and reversible mutations and selection using a Bayesian inference. 
 

## Version 

0.1 (March 2021)


## Citation

In preparation


## Downloading and compiling BEST

First of all, we need to download all the necessary files to compile and run BEST; these files are in this GitHub repository. You can download these files manually or directly from the terminal:


```
git clone https://github.com/mrborges23/best.git best
```

BEST is implemented in C++ and needs to be compiled; this step will produce an executable that can then be used to run BEST. We use the `g++` compiler, but others could be used instead:

```
g++ best.cpp -o best -O3 -std=c++20
```

This step should have created an executable called `best` in your working directory.


## Control file

The control file includes all the necessary parameters to run BEST:
The control file is where you indicate the location of the data, the model to use, and specify the priors to use during Bayesian inference.


| Parameter | Description |
|---|---|
| ```model``` | ```S``` stands for the Moran model with reversible mutations and selection whereas ```N``` stands for the neutral Moran model. |
| ```site_frequency_spectrum_file``` | The location of the sampled site frequency spectrum. |
|```mutation_rate```| The mutation rate. We advise to use mutation rate estimated independently of previous estimate of the effective population size (e.g., pedigree-based, etc.). The problem of using such mutation rates is that the inferences will be circular. Note that Best jointly estimates the effective population size. |
|```population_size_prior_mean```| The prior mean of a normal distribtuion. A good option for a uninformative gamma priors is setting the mean to ... and the variance to ... . The effective population size can only take integer values; we renormalize the normal distribution so that only positive integers are possible to be observed. |
|```population_size_prior_variance```| The prior variance of a normal distribtuion. |
|```fitness_coefficient_prior_mean```| The prior mean of the prior normal distribution. The fiftness coefficients can only take positive values; we renomalized the normal distribution so that it  are possible to be observed. |
|```fitness_coefficient_prior_variance``` | The prior variance of the prior normal distribtuion. |
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


## Running BEST

To run BEST, simply place your sampled site frequecy spectrum and the control file in the folder where the BEST executable is. Then, open the terminal and run the executable `best` followed by the name of the control file:

```
./best control.cf
```

BEST imediatly prints out a short description of the data file. Confirm that this information conforms to your data. To make sure BEST is running, you should get the message `BEST has started!`. When the analyses terminate, you should get the message `BEST has finished!`. BEST periodically writes of the `output_file`.


## Output file

The output file has information about 

## Compilation errors and how to solve them

I hope not!

## Questions and bug reporting

Please use **Issues** to report possible bugs, suggest enhancement features, or if you need help using BEST. If you have more theoretical or biological questions, you can directly contact Rui Borges (ruiborges23@gmail.com).


## License

This program is free software. You can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software. See the GNU General Public License (http://www.gnu.org/licenses/) for more details.
