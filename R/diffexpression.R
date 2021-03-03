#' Make difference expression in cluster with two sample
#' 
#' Make difference bewteen two samples
#' Difference method is \code{\link[Seurat]{FindMarkers}}
#'
#' @param object Seurat object
#' @param design Difference expression, expriment compare with contral, or control comapre experiment
#' @param sample In meta.data of Seurat object, contain your design name.
#' @param min.cells In each cluster, at leat cells in samples. default 10
#' @param ... Other param, see in \code{\link[Seurat]{FindMarkers}}
#' @importFrom Seurat FindMarkers WhichCells
#' @export
MakeDiffSample <- function(object, 
                           design, 
                           sample,
                           min.cells = 10,
                           ...) {
  meta.data <- slot(object = object, name = "meta.data")
  if (!(sample %in% colnames(meta.data))) {
    stop("sample is not in names with object meta.data")
  }
  if (!all(design %in% unique(meta.data[, sample, drop = TRUE]))) {
    stop("Check desing in your sample.")
  }
  cluster <- levels(object)
  diff <- lapply(cluster, FUN = function(x, ...) {
    message(paste0("Find different expression in cluster ", x))
    cells <- WhichCells(object = object, idents = x)
    cells.1 <- cells[cells %in% rownames(meta.data[meta.data[, sample, drop = TRUE] == design[1],])]
    cells.2 <- cells[cells %in% rownames(meta.data[meta.data[, sample, drop = TRUE] == design[2],])]
    if (length(cells.1) < min.cells) {
      message(paste0("In celltype ", x, " sample ", design[1], " cells less ", min.cells))
      return(NULL)
    }
    if (length(cells.2) < min.cells) {
      message(paste0("In celltype ", x, " sample ", design[2], " cells less ", min.cells))
      return(NULL)
    }
    diss <- FindMarkers(object = object,
                        ident.1 = cells.1,
                        ident.2 = cells.2, 
                        ...)
    diss <- tibble::rownames_to_column(diss, "gene")
    return(diss)
  }, ...)
  names(diff) <- cluster
  diff <- diff[!vapply(diff, is.null, FUN.VALUE = logical(1))]
}

#' Annotation cell type
#' 
#' Using reference celltype to annotaton cell
#' 
#' @param object An seurat object
#' @param reference The reference dataset, with rows bing cells, and colums
#' being genes. The last column should be "label"
#' @param k Number of genes to select according to the total entropy difference 
#' cell types, default k = 1000.
#' @return A data.frame, two column.
#'
#' @importFrom scibet SciBet
#' @importFrom Seurat VariableFeatures
#' @export
#load("./database/p11_mouse_brain_cell_label.RData")
#colnames(mouse_p11_barin_count) <- gsub(pattern = " ", replacement = "", colnames(mouse_p11_barin_count))

AnnotationCell <- function(object, reference, k = 1000) {
  if (colnames(reference)[dim(reference)[2]] != "label") 
    stop("The last colname in reference must be label")
  if (!any(rownames(object) %in% colnames(reference)))
    stop("Gene name is not match between reference and object")
  ann <- Matrix::t(object[["RNA"]]@counts[VariableFeatures(mouse), ])
  prd <- scibet::SciBet(reference, ann, k = k)
  cell.type <- data.frame(cell = rownames(ann), predict = prd)
  return(cell.type)
}
