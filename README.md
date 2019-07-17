

The bacollite package offers a set of tools for aligning MALDI spectrum data with 'theroretical' peak positions that have been calculated from sequence and chemistry data. The alignment is carried out by cross-corelation and yields a correlation score and a lag score, both of which can be used to determine the presence or absence of a peptide in the spectrum. 

### Installation

see `INSTALL.md` in this directory

----
# Quickstart guide

## Input data

To carry out a bacollite analysis, you will need:

- A set of peptides OR a collagen sequence (some of these are provided in the package; see "Source Data" below)
- A set of MALDI samples in text or Bruker format, preferably as triplicate samples

The main options regarding input data involve generating appropriate peptides from sequences.

## Processing 

- Initial processing is carried out by the `ms_align` function, which returns correlation and lag scores for each peptide and each sample.
- The correlation and lag scores must then be further processed for classification

## Classification / Interpretation

A variety of approaches exist to carry out classification

- Simple thresholding of correlation and/or lag scores can be used to determine whether a peptide has a 'hit' in the spectrum, which may be sufficient for some analyses
- Alternatively, more sophisticated analysis has been developed with the package, curently the subject of a paper under review. 













----
# Source Data

## Seqeunce Data
### Bioarch Mammalian Collagen Dataset

A version of this dataset is available within the package as an R object called bioarch_mammal_sequences. The data set currently holds sequences for 

* human
* sheep
* cow
* goat

It is based on the  sheet
[https://docs.google.com/spreadsheets/d/1QMIJWZtAZ8zJ4uwbAUpwTp_deBwKeVJ8FoBQopXL_dc/edit?usp=sharing], as at 09:48 on 11 January 2017

### Loading other seqeunce data:

Other sequence sets can be loaded via googleshets: 

```
require(googlesheets)
mc_meta <-gs_title("Mammalian Collagen Sequences v0.0.2")
sequences <- gs_read(mc_meta)
```
Alternatively, the google sheet can be exported as a .csv file to your computer, and then loaded like so: 

```
 bms0.1.2 <- read.table(file="Mammalian Collagen Sequences v0.1.2 - Mammal.csv",sep=",")

```

See [https://docs.google.com/spreadsheets/d/1QMIJWZtAZ8zJ4uwbAUpwTp_deBwKeVJ8FoBQopXL_dc/edit#gid=269004764] for the required format of these sheets



# Getting help with R

Sometimes getting help in R is confusing. Here's some things you can try:

```
bacollite::    # lists availble functions in a pop-out in rstudio
help(foo)      # help about function foo
?foo           # same thing 
apropos("foo") # list all functions containing string foo
example(foo)   # show an example of function foo
```

