#' The 'suez' package.
#' 
#' @description Extends the PANAMA (e)QTL mapping approach to handle repeat measurements and multiple conditions
#' 
#' @docType package
#' @name suez-package
#' @aliases suez
#' @useDynLib suez, .registration = TRUE
#' @import methods
#' @import Rcpp
#' @import rstantools
#' @import doMC
#' @importFrom rstan optimizing
#' 
#' @references 
#' * Fusi N, Stegle O, Lawrence ND. 2012. Joint modelling of confounding factors and prominent genetic regulators provides increased accuracy in genetical genomics studies. PLoS Computational Biology 8:e1002330.
#' 
#' * Knowles* DA, Burrows* CK, Blischak JD, Patterson KM, Serie DJ, Norton N, Ober C, Pritchard JK and Gilad Y (2018), "Determining the genetic basis of anthracycline-cardiotoxicity by molecular response QTL mapping in induced cardiomyocytes", eLife
#' DOI: https://doi.org/10.1371/journal.pcbi.1002330, PMID: 22241974
#' 
#' * Stan Development Team (2018). RStan: the R interface to Stan. R package version 2.17.3. http://mc-stan.org
#' 
NULL
