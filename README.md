PredictDiag
================
predict diagnostic product ions for Hb variants
## Installation

Install from GitHub:

``` r

```
## package

install.packages("seqinr") 
install.packages('OrgMassSpecR') 
library(seqinr) 
library(dplyr) 
library(OrgMassSpecR) 
library(tidyr)

## Steps to get the results
# Input the sequennces of HbA beta and Hb variants (including M at N_termius)
Hbvariants <- read.fasta(file = "Hbvariants.fasta", seqtype = "AA",as.string = FALSE) HbAB <- read.fasta(file = "HbA.fasta", seqtype = "AA",as.string = FALSE)

# Input the possible diagnostic ion list for each AA
diag <- read.csv("finddiag.csv")

# Input reference dignostic ion list of HbA beta
HbA.B <- read.csv("ref mass list_pro_1.csv") %>% select(-c(Ref_Mass, Ref_rel_abun, Variant))

# Run the function
PredictDiag <- PredictDiag(Hbvarinats)
