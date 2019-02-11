# RNAmodR.AlkAnilineSeq
# 
`RNAmodR.AlkAnilineSeq` uses classes and workflows from `RNAmodR` to detect 
`m7G`, `m3C` and `D` in RNASeq data specifically treated according to the 
AlkAnilineSeq protocol.

The vignette from the `RNAmodR.AlkAnilineSeq` package provides details on the
specific analysis of AlkAnilineSeq data. Have a look at the vignette from the
`RNAmodR` package for a general description on the provided workflows.

# Installation

The current version of the RNAmodR.Data package is available from GitHub.

```
remotes::install_github("FelixErnst/RNAmodR.Data")
remotes::install_github("FelixErnst/RNAmodR")
remotes::install_github("FelixErnst/RNAmodR.AlkAnilineSeq")
#
library(RNAmodR.AlkAnilineSeq)
```

A submission to Bioconductor is planned.

# Literature

- Marchand V, Ayadi L, __Ernst FGM__, Hertler J, Bourguignon-Igel V,
Galvanin A, Kotter A, Helm M, __Lafontaine DLJ__, Motorin Y (2018):
"AlkAniline-Seq: Profiling of m7 G and m3 C RNA Modifications at Single
Nucleotide Resolution." Angewandte Chemie (International ed. in English) 57
(51), P. 16785–16790. DOI:
[10.1002/anie.201810946](https://doi.org/10.1002/anie.201810946).