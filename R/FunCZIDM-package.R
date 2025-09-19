#' FunCZIDM: Functional Concurrent Zero-Inflated Dirichlet-Multinomial Regression
#'
#' Implements MCMC for the FunC-ZIDM model and helpers for simulation,
#' fitting, and visualization.
#'
#' ## Getting started
#' - Install: see `README` (GitHub link below).
#' - Main function: `?FunCZIDM`. 
#' - Vignette: see `vignette("FunCZIDMIntro")` (optional, requires `ggtern`), or view online at https://brodyee.github.io/
#' - Full index: `help(package = "FunCZIDM")`
#'
#' ## Links
#' - GitHub: https://github.com/brodyee/FunCZIDM
#'
#' @keywords internal
#' @docType package
#' @name FunCZIDM-package
NULL

## usethis namespace: start
#' @importFrom Rcpp sourceCpp
#' @useDynLib FunCZIDM, .registration = TRUE
## usethis namespace: end
NULL
