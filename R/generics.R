#' Combine multiple data
#'
#' @param object an object
#' @param ... Arguments passed to other methods
#'
#' @rdname merged
#' @export merged
#'
merged <- function(object, ...) {
  UseMethod(generic = "merged", object = object)
}


#' Cells in the quality control
#'
#' @param object An object
#' @param ... Arguments passed to other methods
#'
#' @rdname CellQc
#' @export CellQc
#'
# CellQc <- function(object, ...) {
#   UseMethod(generic = "CellQc", object = object)
# }

setGeneric("CellQc", function(object, ...) object)


#' Running analysis pipeline
#'
#' @param object An object
#' @param variable.features.n Use this many features as variable features after ranking by residual variance; default is 3000
#' @param return.only.var.genes If set to TRUE the scale.data matrices in output assay are subset to contain only the variable genes; default is TRUE
#' @param dims Which dimensions to use as input features, used only if features is NULL; defalut is 1:10
#' @param verbose Controls verbosity; default is TRUE
#' @param ... Arguments passed to other methods
#'
#'
#' @rdname Running
#' @export Running
setGeneric(name = "Running", def = function(object,
                                            variable.features.n = 3000,
                                            return.only.var.genes = TRUE,
                                            dims = 1:10,
                                            resolution = 0.8,
                                            verbose = TRUE,
                                            ...) object )
