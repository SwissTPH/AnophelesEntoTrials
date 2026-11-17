# Outcome probabilities of the basic EHT model

Probabilities of the four hut trial outcomes for one arm, using the same
Markov chain formulas as the basic Stan model (`eht_basic.stan`).

## Usage

``` r
eht_outcome_probs(alpha, mu, pC)
```

## Arguments

- alpha:

  Numeric. Attack rate of the arm.

- mu:

  Numeric. Mortality rate of the arm.

- pC:

  Numeric. Probability of surviving after blood-feeding.

## Value

Named numeric vector of probabilities for UA, UD, FA and FD.
