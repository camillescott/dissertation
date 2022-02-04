# Dissertation

PhD dissertation for UC Davis written with [bookdown](https://bookdown.org/) using a modified [aggiedown](https://github.com/ryanpeek/aggiedown) and [manubot](https://manubot.org/).

The built PDF can be found by:
* Navigating to [Actions](https://github.com/camillescott/dissertation/actions)
* Select most recent commit
* Scroll down to artifacts.

## Local Build

Create the base conda environment:

    mamba env create -f environment.yml
    conda activate dissertation
  
Install the `aggiedown` in the `R` rule environment:

    snakemake -j 1 --use-conda --conda-frontend mamba install_deps
  
Build:

    snakemake -j 1 --use-conda --conda-frontend mamba

The resulting PDF will be in `dissertation.pdf`. The first run will take much longer than subsequent builds while the conda environment is created and `manubot` builds its reference cache.
