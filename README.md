# Functions-for-Staggered-DID-Setting-R-
This fold includes several R functions that can be used for **robustness checks** in **staggered Difference-in-Differences (DID) settings**.
## placebo_treatment
The idea is to randomly assign placebo treatment to units that were never treated in the original dataset. For each treatment year, the number of placebo-treated units is matched to the number of actually treated units in that year. This process is repeated `n` times to generate a distribution of placebo estimates.

By comparing the actual estimated effect to the distribution of placebo estimates, one can evaluate whether the observed treatment effect might be driven by spurious identification or chance patterns.

The document includes two main functions:

- A **TWFE-based** placebo test function using traditional two-way fixed effects estimation.
- A **did_multiplegt-based** placebo test function using the method of de Chaisemartin & D'Haultfœuille for staggered adoption.
