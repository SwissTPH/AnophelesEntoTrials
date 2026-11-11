# AnophelesEntoTrials

AnophelesEntoTrials estimates the entomological efficacy of malaria vector control interventions from entomological trials. It fits Bayesian models of the mosquito feeding cycle to the data with Stan, to estimate intervention effects such as the reduction in host availability, pre-prandial killing and post-prandial killing.

The package provides two predefined models:

| Model | `model_type` | Experimental design | Reference |
|-------|--------------|---------------------|-----------|
| Basic model | `"basic"` | Experimental hut trials (EHT) | Champagne et al. (2025) |
| Semifield model | `"semifield"` | Matched semifield experiments with a spatial repellent | Denz et al. (2021) |

Custom Stan models can also be fitted with the same workflow.

## Installation

```r
install.packages("devtools")

# Install the package from GitHub
devtools::install_github("SwissTPH/AnophelesEntoTrials", build_vignettes = TRUE)
```

The models are compiled with rstan, which requires a C++ toolchain (Rtools on Windows).

## Getting started

```r
library(AnophelesEntoTrials)

# Simulated hut trial data (or import your own with import_ento_data())
eht_data <- format_eht_data(simulate_eht_data(seed = 42))

fit <- fit_ento_model(eht_data, model_type = "basic", seed = 123)
extract_stan_summary(fit, data = eht_data)
```

The vignettes describe the workflow in detail:

- `vignette("format-data")`: importing and formatting entomological data
- `vignette("workflow-predefined-models")`: the predefined models, their mathematical description, fitting, results and posterior predictive checks
- `vignette("workflow-custom-models")`: fitting your own Stan models

## References

- Champagne C, Lemant J, et al. (2025). Cascades of effectiveness of next-generation insecticide-treated nets against malaria, from entomological trials to real-life conditions. *Nature Communications* 16:11162. https://doi.org/10.1038/s41467-025-66130-y
- Denz A, Njoroge MM, et al. (2021). Predicting the impact of outdoor vector control interventions on malaria transmission intensity from semi-field studies. *Parasites & Vectors* 14:64. https://doi.org/10.1186/s13071-020-04560-x
