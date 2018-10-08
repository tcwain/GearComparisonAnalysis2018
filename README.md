# GearComparisonAnalysis2018

This R package contains all data and R-language scripts needed to reproduce the 
analysis in Wainwright et al. (2018, _Effect of a mammal excluder device on 
catches of small pelagic fishes_, **Marine and Coastal Fisheries** _in review_).

The full analysis and detailed results for the publication are in the vignette
"MED_Gear_Analysis.pdf", which can be regenerated on package installation. Data and 
R scripts used in the vignette can be found in the package source subdirectories 
_inst/extdata_ and _inst/scripts_, respectively.

To install the package with an updated vignette in R, run

```
install.packages("devtools")
devtools::install_github("tcwain/GearComparisonAnalysis2018", build_vignettes=TRUE)
```

Building the vignette may require additional system software, including _TeX_ and 
_pandoc_.
