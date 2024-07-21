<!-- README.md is generated from README.Rmd using devtools::build_readme(). Please edit that file -->

[![](https://img.shields.io/badge/lifecycle-stable-green.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)
[![Project Status: Active - The project has reached a stable, usable
state and is being actively
developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![License:
MIT](https://img.shields.io/badge/license-MIT-blue.svg)](https://cran.r-project.org/web/licenses/MIT)
[![](https://img.shields.io/github/last-commit/dputhier/STarlight.svg)](https://github.com/dputhier/STarlight/commits/main)
[![](https://codecov.io/gh/dputhier/STarlight/branch/main/graph/badge.svg)](https://codecov.io/gh/dputhier/STarlight)

# STarlight

STarlight leverages a grid-based strategy via a uniform binning of x/y molecular coordinates for cell segmentation-free analysis of imaging-based spatial transcriptomics data (*e.g.* Merscope, Xenium or CosMx)

<img src="https://github.com/dputhier/STarlight/assets/49205456/6b4da4b3-5be7-40de-95e3-80b851cce8db" style="float: center" class="logo" width="150">

## Installation

### R

The STarlight library is currently not available in CRAN or Bioconductor. To install it from github, use:

    devtools::install_github("dputhier/STarlight")
    library(STarlight)

### Terminal

Download the *tar.gz* from github or clone the main branch. Uncompress and run the following command within STarlight folder:

    R CMD INSTALL .

Then load the library within R.

    library(STarlight)

## Documentation

Documentation is available at
<https://dputhier.github.io/STarlight/>.

## Case study material
Scripts used to produce the Merscope hepatocarcinoma case study figures are available at <https://github.com/dputhier/STarlight_article>. Merscope data used is from [Magen et al.](https://www.nature.com/articles/s41591-023-02345-0).

