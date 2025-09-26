# FunCZIDM Package

This package contains the sampler associated with the "A Bayesian Functional 
Concurrent Zero-Inflated Dirichlet-Multinomial Regression Model with Application
to Infant Microbiome" paper. The package also contains the code to reproduce the
simulations and case study, along with a vignette that walks through how to use 
the package. Lastly, there is a `shiny` app downloadable that can be used to 
easily generate plots.

The layout of this package is as follows:

- `R` and `src` are folders that contain the source code.
- `doc` contains the reference manual and knitted vignette.
- `vignettes` contains the source `.Rmd` file.
- `data` contains the infant covariates and count in an `.rda` object.
- `inst` contains four folders:
  - `extdata` contains the samples obtained from the case study (after burn-in), 
     and the original infant microbiome data from 
     [PNAS](https://www.pnas.org/doi/full/10.1073/pnas.1409497111#data-availability).
  - `scripts` contains the `.R` and `.cpp` scripts to run the case study and 
    simulations.
  - `slurm` contains the slurm scripts to run the simulations and case study on 
    an HPC cluster.
  - `notebooks` contains the .ipynb files to preprocess the data and obtain the
    simulation results.
  - `app` contains the code for the visualization shiny app.

All other files/folders were made automatically and have minimal additions by 
us.

## Vignette (User Guide)

The "Getting Started with FunCZIDM" vignette is hosted 
[here](https://brodyee.github.io/packageVignette/FunCZIDMIntro.html). One can
build it locally; however, `ggtern` is needed to build all the plots. Note, 
`ggtern` is not required to build, as these sections will not evaluate if not 
installed.

## Installation

You can install the development version from GitHub:

```
# install.packages("remotes")
remotes::install_github("brodyee/FunCZIDM", build_vignettes = FALSE)
# Optional: for the vignette and ternary plots
# install.packages("ggtern")
# remotes::install_github("brodyee/FunCZIDM", build_vignettes = TRUE)
```

## Reproducing the Results

We prioritize minimal steps to reproduce our results. We strongly recommend 
reproducing the simulation results with an HPC cluster, as it is very 
computationally intensive. The case study can be run locally, though running it 
on a HPC cluster is preferred. 

If you have access to an HPC cluster with SLURM, the steps to reproduce are as 
follows:

1) Run `Rscript -e "cat(system.file('slurm', package = 'FunCZIDM'))"` to recover
   the path to the `slurm` folder, and copy all the scripts into the folder where 
   you want your results.
2) Fill in the user information for `slurm` in the `.sh` files: email="", ...
3) Run the `.sh`, for example:
   `./simulations.sh`

This will automatically run the `slurm` arrays for the simulation or case study.

If the user does not have access to `slurm`, one can run multiple processes 
simultaneously using the following `bash` script.  

```
#!/usr/bin/env bash
set -euo pipefail

# Get script paths
CASE=$(Rscript -e "cat(system.file('scripts', 'caseStudy.R',
       package='FunCZIDM'))")
TRACE=$(Rscript -e "cat(system.file('scripts', 'chainsTraceplots.R',
        package='FunCZIDM'))")
COMBINE=$(Rscript -e "cat(system.file('scripts', 'combineOutputs.R',
          package='FunCZIDM'))")

# Run four seeds in parallel
Rscript "$CASE" --seed 18342 -i 85000 -b 75000 --toReturn RA \
  --outputFolder "caseStudy" --thin 10 --df 4 --kappaShape 100 --kappaRate 900 &
Rscript "$CASE" --seed 33261 -i 85000 -b 75000 --toReturn RA \
  --outputFolder "caseStudy" --thin 10 --df 4 --kappaShape 100 --kappaRate 900 &
Rscript "$CASE" --seed 87808 -i 85000 -b 75000 --toReturn RA \
  --outputFolder "caseStudy" --thin 10 --df 4 --kappaShape 100 --kappaRate 900 &
Rscript "$CASE" --seed 67235 -i 85000 -b 75000 --toReturn RA \
  --outputFolder "caseStudy" --thin 10 --df 4 --kappaShape 100 --kappaRate 900 &

wait

# Post-process
Rscript "$TRACE"
Rscript "$COMBINE"
```

If using a Windows machine, either using WSL or Git Bash will allow for this 
script to run. There are other options for running multiple R processes 
simultaneously, though they may result in differing results due to the random
number generation.

## `conda` Environment for `shiny` App

While it is not required to make an `conda` environment, we recommend it to 
reduce multiple installs when running directly from `reticulate`. Also, errors
with package installs are easier to sort out within an environment than from 
`reticulate`. Below are the steps to set up the environment.

#### Installing `conda`

If you do not currently have conda installed, you can install it locally 
[here](https://www.anaconda.com/download/success). Choose the appropriate file 
under Miniconda Installers, and click on it after installing.

#### Creating the Environment for the `shiny` app

Everything necessary for the `shiny` app python enviorment is contained in
`inst/app/environment.yml`. One can download it from GitHub or find it locally 
after installing the package with the following line in R:

```
system.file('app', 'environment.yml', package='FunCZIDM')
```

Once you move the file to the working directory, you can create the conda 
environment with the following:

```
conda env create -f environment.yml
```

This will create the `shinyPy` environment used in the shiny app. The app 
will automatically activate the environment, so there is no need to activate it 
after creating it.