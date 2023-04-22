#' @title survMISA : Mediation Identification by splitting and aggregation in survival data
#'
#' @description
#' It is a package designed to implement *survMISA* procedure in practice.
#' *survMISA* procedure is an extension of MISA method propose by Guo et al.(2023).
#' It is a novel procedure, which makes inference without p-values and achieves false discovery rate (FDR) control.
#' By sample splitting and aggregation, we construct statistics with symmetry and deploy the symmetry to obtain a data-driven threshold.
#' We also give a de-randomization method to reduce the effect of sample splitting on the stability of the results.
#'
#' @section About the function:
#' The package supports step-by-step execution of the *survMISA* program, as well as an `AllinOne` function to complete all steps.
#' We also integrate the derandomization operation into `AllinOne` function.
#' To improve operation efficiency, we use `parallel` and `doParallel` package to take parallel computation.
#' This is not necessary, of course, but when the mediation variable dimension is high,
#' doing all the steps sequentially consumes a lot of computation time.
#'
#' @references
#' Guo, X., Ren, H., Zou, C., and Li, R. (2023). Large-scale mediation effect signal detection. Working paper.
#'
#' --survMISA--
#'
#' @docType package
#' @name survMISA

NULL


