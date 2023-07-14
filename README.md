# microbiome-inference

This repository includes the code used to produce the results presented 
in *Zapién-Campos, et al.*, 2023.

The aim is to use Bayesian inference to estimate microbiome data interactions
using theoretical stochastic models,  extracting the information contained 
in the statistical moments of microbial abundance.

## Organization

1. Set Environment (**microbiome-inference/install_environment**).

Code to install the Anaconda environment used by the authors.
All code is written in Python 3.8.

2. Data (**microbiome-inference/data**).

Directory storing empirical data and to store simulated data.

3. Figures (**microbiome-inference/figures**).

Scripts to produce figures 2 and 3.

4. Simulations and Numerics (**microbiome-inference/simulations_n_numerics**).

Scripts to simulate microbiomes using the Gillespie algorithm and 
to solve models numerically.

## Usage

Code files are structured in the following way:

`sc_*.py` contains the source code.

`par_*.py` contains the parameters.

`exe_*.py` produces the output.

`*.ipynb` IPython notebooks to generate and analyse data.
Prefix numbers indicate the suggested order of execution 
of the individual notebooks.