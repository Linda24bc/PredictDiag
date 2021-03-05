PredictDiag
================
predict diagnostic product ions for Hb variants
## Installation

Install from GitHub:

``` 

```
## Introduction
+ This procedure can predict the diagnostic ions of Hb beta variants that have been experimentally teseted by refering to the product ions of Hb A beta.
+ The output list from step 5 can be combined with the original databse in Variants Identifier as the updated database for sereaching.

### Prerequisite 1: Install packages
+ install.packages("seqinr") 
+ library(seqinr) 
+ library(dplyr) 
+ library(tidyr)
+ library(tidyverse)
+ library(data.table)

### Prerequisite 2.1: Load monomz and PredictDiag function
The R apckage is avaiable in Github. 

monomz(sequence, fragments = "by") # for cz ions, define fragments as "cz"

PredictDiag(Hbvarinats)

### Prerequisite 2.2: Input the list of residue numbers of possible diagnostic ions for each AA in the Hb beta

diag <- read_csv("finddiag.csv")

### Prerequisite 2.3: Input list of reference product ions for HbA beta

HbA.B <- read_csv("ref mass list_pro_1.csv") %>% select(-c(Ref_Mass, Ref_rel_abun, Variant))

### Step 1: Input the sequennces of HbA beta and Hb beta variants 
Mutiple sequences of variants sequences can be included in one .fasta file, the sequences should have the N-terminal Met while the comparision results exclude the N-ternimal Met.

Hbvariants <- read.fasta(file = "Hbvariants_new.fasta", seqtype = "AA",as.string = FALSE)

HbAB <- read.fasta(file = "HbA.fasta", seqtype = "AA",as.string = FALSE)

### Step 2: Predict the diagnostic ions by running the function

PD.result <- PredictDiag(Hbvarinats)

### Step 3:Output results in .csv file

write.csv(PD.result, "PredictDiag_new.csv", row.names = FALSE)
