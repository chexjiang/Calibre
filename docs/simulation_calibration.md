# Simulation calibration findings

Target FDR is 0.10. Formal Calibre simulations use the standardized Pearson
residual bootstrap (`null_generation = "s"`) together with both finite-sample
components:

1. the multiplicative correction
   `1 / (1 + sqrt(log2(n_pairs) / n_cells))`; and
2. the plus-one numerator in the empirical FDP estimator.

The plus-one term is necessary. With the earlier `q / 2` numerator, an isolated
largest observed statistic is accepted whenever no synthetic statistic exceeds
it, which occurs with probability near one half under the global null. A
multiplicative adjustment to `q` cannot fix that event.

## High-replicate calibration diagnostic

The `diagnostic` profile uses 500 pairs, 200 cells, and 100 replicates for each
Poisson, NB, and ZINB DGP at effect sizes 0 and 0.2. All 3,000 method fits were
valid. Under the global null, corrected Calibre-Poisson and Calibre-NB made no
discoveries in any of the 600 method-specific replicate runs. The uncorrected
archived run had discovery probabilities between 0.44 and 0.55.

## Weak-signal power pilot

The `power_pilot` profile uses 500 pairs, 500 cells, a 20% signal fraction,
effect size 0.2, and 20 replicates per DGP. Corrected Calibre controlled FDR in
all six DGP/model combinations, but it did not have uniformly greater power
than scMultiMap. This pilot must not be presented as evidence of absolute
superiority.

Constant, ATAC-independent zero inflation does not invalidate the NB mean
coefficient: `E(Y | ATAC) = (1 - pi) * mu`, so a constant `pi` changes only the
intercept on the log scale. Consequently, NB-GLM need not lose FDR control under
the current ZINB DGP. An additional ZINB stress test should be motivated by a
specific biological dropout mechanism and must preserve a valid null rather
than making dropout directly ATAC-dependent for pairs labelled null.

## All-method effect-size trend

The `trend` profile adds Signac and evaluates ten effect sizes from 0 to 1.2
with 500 cells, 500 candidate pairs, a 20% signal fraction, and 20 replicates
per DGP/effect combination. It contains 600 simulated datasets and 3,600
method-level runs; every method returned 500 valid candidate-pair results in
every dataset.

The Signac benchmark follows the preprocessing and testing choices in
`incoming/*/run_Signac.R`: RNA `NormalizeData`, ATAC counts, Pearson
association, matched background peaks based on GC content/count/sequence
length, two-sided p-values, and BH over the common candidate-pair family. The
project contains both a direct `Signac::LinkPeaks` adapter and a
candidate-restricted implementation of the same statistic. The latter is used
for the high-replicate trend because it avoids repeatedly testing incidental
cis pairs outside the common benchmark universe.

Across nonzero effects, the two corrected Calibre methods had mean FDR values
between approximately 0.049 and 0.063 after averaging within each DGP/effect
combination. Signac power was lower for weak effects, particularly under ZINB,
and approached one as the effect increased. Global-null estimates based on 20
replicates remain coarse and should be interpreted together with the separate
100-replicate calibration profile.

## Sample-size trend

The `model_based` profile repeats the all-method trend at 100, 200, 300, and
400 cells, using exactly the same 500 pairs, 20% signal fraction, ten effect
sizes, three DGPs, 20 replicates, and six methods as the 500-cell trend. The new
run contains 2,400 datasets and 14,400 method-level rows with no method errors.
It is combined with the existing 500-cell run only during summarization.

Corrected Calibre made no global-null discoveries at any tested cell number.
Across nonzero effects, mean FDR ranges over DGP/effect combinations remained
below 0.068 for both Calibre versions. Power increased strongly with cell
number, especially near the detection boundary. For example, at effect 0.4,
Calibre-NB power under the NB DGP increased from 0.037 at 100 cells to 0.941 at
500 cells; under ZINB it increased from 0 to 0.554. Calibre-Poisson under ZINB
increased from 0.042 to 0.694 at the same effect.

Two 100-cell datasets contained one peak that did not pass the existing
Calibre/GLM variability filter, giving 499 rather than 500 fitted pairs for the
four regression-based methods. All other regression runs and every scMultiMap
and Signac run fitted all 500 pairs.
