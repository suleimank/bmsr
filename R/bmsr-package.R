#' The 'bmsr' package.
#'
#' @description The bmsr package implements joint regression from multiple data sources in a Bayesian framework. The package provides implementation for both single-task and multi-task regression. The model is implemented using STAN and interface is provided using R programming language. Options for training the model using both NUTS sampler and variational inference are provided. The package is structured for ease of use and the included demo shows the model execution on real-life as well as simulated datasets.
#'
#' @docType package
#' @name bmsr-package
#' @aliases bmsr
#' @useDynLib bmsr, .registration = TRUE
#' @import methods
#' @import Rcpp
#' @importFrom rstan sampling
#'
#' @references
#' Brian S. White, Suleiman A. Khan, Mike J Mason, Muhammad Ammad-ud-din, Swapnil Potdar, Disha Malani, Heikki Kuusanmäki, Brian J. Druker, Caroline A Heckman, Olli Kallioniemi, Stephen E Kurtz, Kimmo Porkka, Cristina E. Tognon, Jeffrey W. Tyner, Tero Aittokallio, Krister Wennerberg, Justin Guinney,
#' \emph{Bayesian multi-source regression and monocyte-associated gene expression predict BCL-2 inhibitor resistance in acute myeloid leukemia},
#' <To Appear>, (2021)
NULL
