# Functions for Staggered DID Setting (R)
This fold includes several R functions that can be used for **robustness checks** in **staggered Difference-in-Differences (DID) settings**.
## placebo_treatment
The idea is to randomly assign placebo treatment to units that were never treated in the original dataset. For each treatment year, the number of placebo-treated units is matched to the number of actually treated units in that year. This process is repeated `n` times to generate a distribution of placebo estimates.

By comparing the actual estimated effect to the distribution of placebo estimates, one can evaluate whether the observed treatment effect might be driven by spurious identification or chance patterns.

The document includes two functions:

- A **TWFE-based** placebo test function using traditional two-way fixed effects estimation.
- A **did_multiplegt-based** placebo test function using the method of de Chaisemartin & D'Haultfœuille for staggered adoption.

## match_single_variable
This function constructs matched control groups for staggered treatment adoption based on nearest-neighbor matching using the Euclidean distance of a single outcome variable. 

For each treatment year T, treated units are matched to never-treated units using pre-treatment outcomes from $T−k$ to $T−m$ (user-specified). Matching quality is evaluated using the average absolute difference in outcomes during $T−q$ to $T−1$ (user-specified), and matches are retained only if this difference is below a user-defined tolerance.

The output is a list of matched datasets for each treatment cohort.
