# microbiome-inference

This repository includes the code used to produce the results presented 
in `Zapién-Campos R, et al., 2024`.

The aim is to use Bayesian inference to estimate microbiome data interactions
using stochastic modelling,  extracting the information contained 
in the statistical moments of microbial abundance.

## Organization

1. Set Environment (**microbiome-inference/install_environment/**).

Code to install the Anaconda environment used by the authors.
pyABC, used to perform Approximate Bayesian Computation (ABC) must be installed within the created Anaconda environment as indicated in https://github.com/romanzapien/microbiome-inference.
All Python code is written in Python 3.8.

2. Structural Identifiability (**microbiome-inference/identifiability/**)

Python scripts to make the Matlab files to perform structural identifiability analyses using GenSSI 2.0 (https://github.com/genssi-developer/GenSSI).

3. Code to make Figures (**microbiome-inference/fig*/**)

    3.1 Data (**microbiome-inference/fig*/data/**).

    Directory to store simulated or empirical data.

    3.2 Figures (**microbiome-inference/fig*/figures/**).

    Jupyter notebooks to make the figures.

    3.3 Simulations and Numerics (**microbiome-inference/fig*/simulations_n_numerics/**).

    Jupyter notebooks and scripts to simulate using the Gillespie algorithm, to solve ODE models numerically, and to perform ABC inference.

## Usage

The Python code files are structured in the following way:

`sc_*.py` contains the source code.

`par_*.py` contains the parameters.
The content within these files must be commented/uncommented depending on the `exe_*.py` file to execute. Specific instructions are given within the file.

`exe_*.py` produces the output.

`*.ipynb` Jupyter notebooks to generate and analyse data.
Prefix numbers indicate the suggested order of execution 
of the individual notebooks.

The Matlab-related code files are structured in the following way:

`make_*.py` generates the Matlab files.
Within these `make_*.py` files users can select the model (logistic or Lotka-Volterra), dynamics (deterministic or stochastic), and number of microbial types.

`*_abund.m` constains the generated Matlab code for a given model.

`*_abund_run.m` constains the executable Matlab code for a given model.

The output of the structural identifiability analyses are stored in the directories `run*/`.