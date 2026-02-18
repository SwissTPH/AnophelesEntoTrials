# Install the package
```
install.packages("devtools") #Prerequisites

# Install the package from GitHub
devtools::install_github('SwissTPH/AnophelesEntoTrials', build_vignettes = TRUE)
```
## Explore vignettes to start using the package

You can utilize default models defined in the package or use your own one. 

Be careful while importing the data, to match your parameter-specific inputs.

## To-Dos
-> More testing with different data structures to ensure edge-cases handling

-> Test and update semi-field experiments model

-> Split fit_models.R into parameter set-up script and main model-calling script to make it easier to maintain

-> Ensure right parameter init function OR communicate with user more explicitly which parameters to use

-> Update vignettes from unclear moments
