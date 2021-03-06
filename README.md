# RNAmodR.AlkAnilineSeq [![Build Status](https://travis-ci.com/FelixErnst/RNAmodR.AlkAnilineSeq.svg?branch=master)](https://travis-ci.com/FelixErnst/RNAmodR.AlkAnilineSeq) [![codecov](https://codecov.io/gh/FelixErnst/RNAmodR.AlkAnilineSeq/branch/master/graph/badge.svg)](https://codecov.io/gh/FelixErnst/RNAmodR.AlkAnilineSeq)

`RNAmodR.AlkAnilineSeq` uses classes and workflows from `RNAmodR` to detect 
`m7G`, `m3C` and `D` in RNASeq data specifically treated according to the 
AlkAnilineSeq protocol.

The vignette from the `RNAmodR.AlkAnilineSeq` package provides details on the
specific analysis of AlkAnilineSeq data. Have a look at the vignette from the
`RNAmodR` package for a general description on the provided workflows.

# Installation

The current version of the RNAmodR.AlkAnilineSeq package is available from
Bioconductor.

```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("RNAmodR.AlkAnilineSeq")
library(RNAmodR.AlkAnilineSeq)
```

# Literature

- Marchand, Virginie, Lilia Ayadi, __Felix G. M. Ernst__, Jasmin Hertler,
Valérie Bourguignon-Igel, Adeline Galvanin, Annika Kotter, Mark Helm, 
__Denis L. J. Lafontaine__, and Yuri Motorin. 2018. “AlkAniline-Seq: Profiling 
of m7G and m3C Rna Modifications at Single Nucleotide Resolution.” Angewandte 
Chemie International Edition 57 (51): 16785–90. 
https://doi.org/10.1002/anie.201810946.
