# MSstatsTMT

<!-- badges: start -->
[![Bioconductor Release Build](https://bioconductor.org/shields/build/release/bioc/MSstatsTMT.svg)](https://bioconductor.org/checkResults/release/bioc-LATEST/MSstatsTMT/)
[![Codecov test coverage](https://codecov.io/gh/Vitek-Lab/MSstatsTMT/branch/devel/graph/badge.svg)](https://codecov.io/gh/Vitek-Lab/MSstatsTMT/branch/devel)
[![Bioconductor Downloads Rank](https://bioconductor.org/shields/downloads/release/MSstatsTMT.svg)](https://bioconductor.org/packages/stats/bioc/MSstatsTMT/)
[![Years in Bioconductor](https://bioconductor.org/shields/years-in-bioc/MSstatsTMT.svg)](https://bioconductor.org/packages/release/bioc/html/MSstatsTMT.html#since)
[![License: Artistic-2.0](https://img.shields.io/badge/license-Artistic--2.0-blue.svg)](https://opensource.org/licenses/Artistic-2.0)
<!-- badges: end -->

MSstatsTMT is an R/Bioconductor package for statistical relative quantification
of proteins and peptides in shotgun mass spectrometry-based proteomics
experiments with isobaric labeling. It is applicable to tandem mass tag (TMT)
and iTRAQ workflows, and is designed for designs with multiple mixtures and
multiple technical replicate MS runs per mixture, including case-control and
repeated-measures (time-course, paired) designs. Given peptide-spectrum matches
quantified by upstream tools (Proteome Discoverer, MaxQuant, OpenMS, SpectroMine,
FragPipe/Philosopher), MSstatsTMT performs normalization (including
reference-channel normalization), protein-level summarization, and model-based
statistical testing based on linear mixed-effects models to detect
differentially abundant proteins across conditions.

MSstatsTMT is part of the [MSstats](https://github.com/Vitek-Lab/MSstats)
family of packages, developed and maintained by the
[Vitek Lab](https://olga-vitek-lab.khoury.northeastern.edu/) at Northeastern
University. The package and its documentation are also available at
[msstats.org/msstatstmt](http://msstats.org/msstatstmt/).

## Installation

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("MSstatsTMT")
```

The development version can be installed directly from this repository:

```r
BiocManager::install("Vitek-Lab/MSstatsTMT", ref = "devel")
```

## Quick Start

```r
library(MSstatsTMT)

# Example Proteome Discoverer output, already converted to MSstatsTMT format
data(input.pd)

# Protein-level summarization with global and reference-channel normalization
quant <- proteinSummarization(input.pd,
                              method = "msstats",
                              global_norm = TRUE,
                              reference_norm = TRUE)

# Test all pairwise comparisons between conditions
test <- groupComparisonTMT(quant, contrast.matrix = "pairwise", moderated = TRUE)

head(test$ComparisonResult)
```

Real experiments typically start by converting a search tool's output into
MSstatsTMT format with one of the `*toMSstatsTMTFormat()` converters below, then
proceed with `proteinSummarization()` and `groupComparisonTMT()` as above. See
the vignette for a complete, worked example.

## Supported Converters

MSstatsTMT does not read raw search-tool output directly. Instead, a converter
translates each tool's report into MSstatsTMT format before
`proteinSummarization()` is called:

| Search tool / format | Converter function |
| --- | --- |
| Proteome Discoverer | `PDtoMSstatsTMTFormat()` |
| MaxQuant | `MaxQtoMSstatsTMTFormat()` |
| SpectroMine | `SpectroMinetoMSstatsTMTFormat()` |
| OpenMS | `OpenMStoMSstatsTMTFormat()` |
| Philosopher / FragPipe | `PhilosophertoMSstatsTMTFormat()` |

See the [MSstatsTMT vignette](vignettes/MSstatsTMT.Rmd) for the required input
files and options for each converter.

## Documentation

- [MSstatsTMT vignette](vignettes/MSstatsTMT.Rmd) — full worked example from raw data to results
- [Official website: msstats.org/msstatstmt](http://msstats.org/msstatstmt/)
- [Bioconductor package page and reference manual](https://bioconductor.org/packages/MSstatsTMT)

## Getting Help / Reporting Bugs

- **Questions about usage, statistical methods, or troubleshooting:** please
  post to the [MSstats Google Group](https://groups.google.com/forum/#!forum/msstats).
  This is monitored by the development team and searchable, so it's the fastest
  way to get help and to see if your question has already been answered.
- **Bug reports and feature requests for this repository:** please open a
  [GitHub issue](https://github.com/Vitek-Lab/MSstatsTMT/issues).

## References

If you use MSstatsTMT, please cite the relevant publication(s):

1. Huang T, Choi M, Tzouros M, Golling S, Pandya NJ, Banfai B, Dunkley T,
   Vitek O. **MSstatsTMT: Statistical Detection of Differentially Abundant
   Proteins in Experiments with Isobaric Labeling and Multiple Mixtures.**
   *Mol Cell Proteomics*. 2020;19(10):1706-1723.
   [DOI: 10.1074/mcp.RA120.002105](https://doi.org/10.1074/mcp.RA120.002105)
2. Huang T, Staniak M, da Veiga Leprevost F, Figueroa-Navedo AM, Ivanov AR,
   Nesvizhskii AI, Choi M, Vitek O. **Statistical Detection of Differentially
   Abundant Proteins in Experiments with Repeated Measures Designs and Isobaric
   Labeling.** *J Proteome Res*. 2023;22(8):2641-2659.
   [DOI: 10.1021/acs.jproteome.3c00155](https://doi.org/10.1021/acs.jproteome.3c00155)
3. Figueroa-Navedo AM, Kapre R, Gupta T, Xu Y, Phaneuf CG, Jean Beltran PM,
   Xue L, Ivanov AR, Vitek O. **MSstatsTMT Improves Accuracy of Thermal Proteome
   Profiling.** *Mol Cell Proteomics*. 2025;24(8):100999.
   [DOI: 10.1016/j.mcpro.2025.100999](https://doi.org/10.1016/j.mcpro.2025.100999)

## Funding

MSstats development has been supported by the Chan Zuckerberg Initiative's
[Essential Open Source Software for Science](https://chanzuckerberg.com/eoss/proposals/).

## License

MSstatsTMT is released under the [Artistic-2.0](https://opensource.org/licenses/Artistic-2.0) license.
