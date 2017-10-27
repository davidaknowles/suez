# suez: an R package for mapping response-eQTLs

`suez` extends the `PANAMA` eQTL mapping framework (http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1002330) to handle
- multiple conditions (e.g. stimulated/unstimulated), including response-eQTL mapping
- latent factor correction
- known kinship between individuals

## Installation 

The easiest way to install is
```
if (!require("devtools")) install.packages("devtools", repos='http://cran.us.r-project.org')
devtools::install_github("davidaknowles/suez")
```

## Usage 

For example usage see [here](https://github.com/davidaknowles/dox/blob/523fb973e06daad34ce621df6a7644df27be02e2/code/panama_test.R). 
