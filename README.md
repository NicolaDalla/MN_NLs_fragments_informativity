# Informative Neutral Losses and Fragment Ions in Molecular Networks

This repository provides an R implementation of a network-based framework to rank neutral losses and fragment ions according to their relevance for molecular network interpretation in untargeted metabolomics.

Rather than modifying spectral similarity or network construction, this method adds an interpretative layer that quantifies how strongly individual spectral features (neutral losses or fragment m/z values) are associated with specific regions of a molecular network.

The approach is based on network modularity, permutation testing, and frequency-adjusted ranking using Quantile Generalized Additive Models (QGAMs).

---

## Concept

In molecular networking, clusters of related MS/MS spectra can be visualized, but it is often unclear which neutral losses or fragment ions actually define these clusters. Common or ubiquitous features may dominate similarity calculations without providing diagnostic structural information.

This framework identifies spectral features whose presence is non-randomly structured across the network and therefore likely to be chemically informative.

For each spectral feature:
- nodes are split into those containing the feature and those that do not,
- network modularity is computed for this partition,
- a permutation-based Z-score is calculated,
- frequency-dependent bias is removed using QGAM,
- features are ranked by their QGAM residuals.

---

## Main functions

A main function is provided:

### `MN_info()`
Ranks **neutral losses** or **fragment ion m/z values**  by their informativeness in a molecular network.

Function return a data frame where spectral features are ordered by their QGAM residuals, providing a prioritized list of features most characteristic of specific molecular network regions.

### `Real_application.Rmd`

contains the example of application of the function to an Alternaria fungi extracts datasets acquired in DDA.  
feature creation, MS/MS extraction and plotting code of MN are included.

---

## Input

The function require:
- "method" parameter: "FIs" if we want to investigate Fragment Ions informativity, "NLs" if we want to investigate neutral losses informativity.
- an MS/MS dataset in **Spectra** format,
- a similarity matrix or molecular network representation.

The implementation allows users to control:
- DBscan grouping parameters
- spectral similarity threshold for network construction,
- number of permutations,
- m/z tolerance,
- minimum and maximum feature prevalence.

---

## Output

Each function returns a data frame containing, for each neutral loss or fragment ion:
- feature m/z,
- number of spectra containing the feature,
- modularity,
- permutation modularity-based Z-score,
- QGAM residual (informativeness score).

Higher residuals indicate features that are more specifically associated with molecular network substructures.

---

## Dependencies

This package relies on:

- `Spectra` for MS/MS data handling  
- `igraph` for molecular network analysis  
- `dbscan` for m/z grouping  
- `qgam` for Quantile Generalized Additive Models  

---

## Typical workflow

1. Import MS/MS data as a `Spectra` object and calculte spectral similarities 
2. Construct a molecular network  
3. Run `MN_info()` for FIs or NLs "method" parameter (the merging of the two output is possible)
4. Visualize feature informativeness using the returned ranking or built-in plotting functions  
