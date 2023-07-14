#! /bin/sh

unset PYTHONPATH
unset PYTHONHOME

conda env create --file environment.yml
conda activate ABCSMC

ipython kernel install --user --name ABCSMC
