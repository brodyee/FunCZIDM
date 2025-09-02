# FunCZIDM Package

This package contains the sampler associated with the "A Bayesian Functional 
Concurrent Zero-Inflated Dirichlet-Multinomial Regression Model with Application
to Infant Microbiome" paper. The package also contains the code to reproduce the
simulations and case study, along with a vignette that walks through how to use 
the package. Lastly, there is an associated `shiny` app downloadable 
![here](willBeHere).

The layout of this package is as follows:

- `R` and `src` are folders that contain the source code.
- `doc` contains the reference manual and knitted vignette.
- `vignettes` contains the source `.Rmd` file.
- `data` contains the infant covariates and count in an `.rda` object.
- `inst` contains four folders:
  - `extdata` contains the samples obtained from the case study (after burn-in), 
     and the original infant microbiome data from 
     ![PNAS](https://www.pnas.org/doi/full/10.1073/pnas.1409497111#data-availability).
  - `scripts` contains the `.R` and `.cpp` scripts to run the case study and 
    simulations.
  - `slurm` contains the slurm scripts to run the simulations and case study on 
    an HPC cluster.
  - `notebooks` contains the .ipynb files to preprocess the data and obtain the
    simulation results.

All other files/folders were made automatically and have minimal additions by 
us.

### Reproducing the Results

We prioritize minimal steps to reproduce our results. We strongly recommend 
reproducing the simulation results with an HPC cluster, as it is very 
computationally intensive. The case study can be run locally, though running it 
on a HPC cluster is preferred. 

If you have access to an HPC cluster with SLURM, the steps to reproduce are as 
follows:

1) Run `Rscript -e "cat(system.file('slurm', package = 'FunCZIDM'))"` to recover
   the path to the `slurm` folder, and copy all the scripts into the folder you 
   want your results.
2) Fill in the user information for slurm in the `.sh` files: email="", ...
3) Run the `.sh`, for example:
   `./simulations.sh`

This will automatically run the slurm arrays for the simulation or case study.

If you do not have access to an HPC cluster, one can follow the steps in the 
vignette to run the sampler for the case study data. For multiple chains, one 
can run multiple processes simultaneously, use screen in Linux, or run the 
chains sequentially.