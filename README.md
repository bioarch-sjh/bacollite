---
output: html_document
---


#Getting help

##Vignettes



#Seqeunce Data

##Bioarch Mammalian Collagen Dataset

A version of this dataset is available within the package as an R object called bioarch_mammal_sequences. The data set currently holds sequences for 

* human
* sheep
* cow
* goat

Other sequence sets can be loaded via googleshets: 

```
require(googlesheets)
mc_meta <-gs_title("Mammalian Collagen Sequences v0.0.2")
sequences <- gs_read(mc_meta)
```

See [https://docs.google.com/spreadsheets/d/1QMIJWZtAZ8zJ4uwbAUpwTp_deBwKeVJ8FoBQopXL_dc/edit#gid=269004764] for the required format of these sheets



