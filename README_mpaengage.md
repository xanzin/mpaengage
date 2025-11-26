# mpaengage

## Overview
This repository contains data and R code to calculate and explore a marine protected area (MPA) engagement and vulnerability index. The workflow standardises input data, applies a vulnerability model, and provides a Shiny application for interactive inspection of results.

## Repository structure
- `model.R` – Main script to run the vulnerability/engagement model.
- `R/` – Supporting R functions.
- `data/` – Core processed data used by the model.
- `inputs/` – Input tables (e.g. indicators, weights, classifications).
- `outputs/` – Model outputs (vulnerability and engagement metrics).
- `figures/` – Saved plots and figures.
- `shiny/Vulnerability_Index/` – Shiny app for interactive exploration of the index.

## Usage
Clone the repository and set the working directory to the project root:

```bash
git clone https://github.com/xanzin/mpaengage.git
cd mpaengage
```

Run the main model script in R:

```r
source("model.R")
```

This will read the input data, run the index calculations, and write results to `outputs/` (and, where implemented, update figures in `figures/`).

## Shiny application
To launch the Shiny app for interactive exploration of the vulnerability index:

```r
setwd("shiny/Vulnerability_Index")
shiny::runApp()
```

The app uses the model outputs to visualise spatial patterns, indicator scores, and summary statistics for MPAs and associated social–ecological units.

## Requirements
- R (>= 4.0)
- Standard R packages for data manipulation, plotting, and Shiny (e.g. `tidyverse`, `sf`, `shiny`, `shinydashboard`, `ggplot2`).  
Install missing packages as needed, for example:

```r
install.packages(c("tidyverse", "sf", "shiny", "shinydashboard", "ggplot2"))
```

## Notes
The repository is intended as a reproducible framework for MPA engagement and vulnerability assessments. It can be adapted to new case studies by replacing the contents of `inputs/` and, where required, updating the model code in `model.R` and `R/`.

## License
To be specified.
