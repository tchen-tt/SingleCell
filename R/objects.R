#' @include generics.R
#' @importClassesFrom Seurat Seurat
#'
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Class definitions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



#' The DataSheet Class
#'
#' The DataSheet class is an intermediate data storage class that sotres the data path
#'
#' @slot data.path Vector of data path
#' @slot sample.name Vector of data name
#'
#' @name DataSheet-class
#' @rdname DataSheet-class
#' @exportClass DataSheet
#'
DataSheet <- setClass(
  Class = "DataSheet",
  slots = list(
    data.path = "vector",
    sample.name = "vector"
  )
)





#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Create DataSheet
#'
#' Use input data information create a DataSheet object
#'
#' @param data.path Vector of data path
#' @param sample.name Vector of data name
#'
#' @return DataSheet object
#' @export
#'
#' @examples
#' data <- CreateSingleRNA(data.path = c("./data/filtered_gene_bc_matrices3k/hg19/", "./data/filtered_feature_bc_matrix10k/"),
#'                         sample.name = c("pbmc", "lymphoma"))
#'
CreateSingleRNA <- function(data.path, sample.name) {
  if (length(data.path) != length(sample.name)) {
    stop("In put data path is not match to your input sample name.")
  }
  return(new(Class = "DataSheet", data.path = data.path, sample.name = sample.name))
}



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Methods for singleRNA-define generics
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @importFrom Seurat Read10X
#' @importFrom Seurat CreateSeuratObject
#'
#' @rdname merged
#' @export
#' @method merged DataSheet
#'

merged.DataSheet <- function(object, ...) {
  single <- parallel::mclapply(1L:length(object@data.path), FUN = function(x, ...) {
    counts <- Read10X(data.dir = object@data.path[x])
    counts <- CreateSeuratObject(counts = counts, ...)
    counts[['sample']] <- object@sample.name[x]
    return(counts)
  }, ...)
  if (length(single) == 1L) {
    return(single[[1L]])
  } else {
    single <- merge(x = single[[1L]], y = single[[-1L]], add.cell.ids = object@sample.name)
  }
  return(single)
}

#' @param min.nFeature_RNA Include cells where at least this many features are detected, default 200
#' @param max.nFeature_RNA Include cells where limits this many features are detected, defalult Inf
#' @param max.percent.mt Include cells where largest proportion of mitochondrial genes, default five percent
#'
#' @importFrom Seurat PercentageFeatureSet DefaultAssay
#'
#' @return Returns the filtered Seruat object.
#'
#' @rdname CellQc
#' @method  CellQc Seurat
#'
# CellQc.Seurat <- function(object,
#                           min.nFeature_RNA = 200,
#                           max.nFeature_RNA = Inf,
#                           max.percent.mt = 5, ...) {
#   assay <- DefaultAssay(object)
#   features <- rownames(object[[assay]]@counts)
#   if (is.null(grep(pattern = "^MT-|^mt-", features))) {
#     warning("No mitochondrial genes were detected. Check if the gene annotation file contains mitochondrial genes",
#             call. = FALSE, immediate. = TRUE)
#   }
#   object[['percent.mt']] <- PercentageFeatureSet(object = object, pattern = "^MT-|^mt-")
#   object <- subset(object, subset = (nFeature_RNA > min.nFeature_RNA) & (nFeature_RNA < max.nFeature_RNA) &
#                      (percent.mt < max.percent.mt))
# }


setMethod("CellQc", "Seurat", function(object,
                                       min.nFeature_RNA = 200,
                                       max.nFeature_RNA = Inf,
                                       max.percent.mt = 5, ...) {
  assay <- DefaultAssay(object)
  features <- rownames(object[[assay]]@counts)
  if (is.null(grep(pattern = "^MT-|^mt-", features))) {
    warning("No mitochondrial genes were detected. Check if the gene annotation file contains mitochondrial genes",
            call. = FALSE, immediate. = TRUE)
  }
  object[['percent.mt']] <- PercentageFeatureSet(object = object, pattern = "^MT-|^mt-")
  object <- subset(object, subset = (nFeature_RNA > min.nFeature_RNA) & (nFeature_RNA < max.nFeature_RNA) &
                     (percent.mt < max.percent.mt))
})


#' @rdname Running
#' @method Running Seurat
#'
setMethod("Running", signature = "Seurat",
          definition = function(object,
                                variable.features.n = 3000,
                                return.only.var.genes = TRUE,
                                dims = 1:10,
                                resolution = 0.8,
                                verbose = TRUE,
                                ...) {
  object <- SCTransform(object = object,
                        variable.features.n = variable.features.n,
                        return.only.var.genes = return.only.var.genes,
                        verbose = verbose, ...)
  object <- RunPCA(object = object,
                   verbose = verbose,
                   ...)
  object <- RunTSNE(object = object,
                    verbose = verbose,
                    dims = dims,
                    ...)
  object <- RunUMAP(object = object,
                    dims = dims,
                    verbose = verbose,
                    ...)
  object <- FindNeighbors(object = object,
                          verbose = verbose,
                          ...)
  object <- FindClusters(object = object,
                         resolution = resolution,
                         verbose = verbose,
                         ...)
  return(object)
  })








