

# Installing bacollite

The bacollite package was built using devtools, and is hosted on github. To install the package in R, you'll need to have devtools installed first:

```
install.packages("devtools")
```


Then load the library

```
library(devtools)
```

Then install bacollite:

```
devtools::install_github("bioarch-sjh/bacollite")
```

# Linux-Ubuntu notes

If devtools install fails, you may need to install curl, which has to be done outside of R in a Terminal. Type the following: 

```
sudo apt-get install libcurl4-openssl-dev
```
