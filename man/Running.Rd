% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generics.R, R/objects.R
\name{Running}
\alias{Running}
\alias{Running,Seurat-method}
\title{Running analysis pipeline}
\usage{
Running(
  object,
  variable.features.n = 3000,
  return.only.var.genes = TRUE,
  dims = 1:10,
  resolution = 0.8,
  verbose = TRUE,
  ...
)

\S4method{Running}{Seurat}(
  object,
  variable.features.n = 3000,
  return.only.var.genes = TRUE,
  dims = 1:10,
  resolution = 0.8,
  verbose = TRUE,
  ...
)
}
\arguments{
\item{object}{An object}

\item{variable.features.n}{Use this many features as variable features after ranking by residual variance; default is 3000}

\item{return.only.var.genes}{If set to TRUE the scale.data matrices in output assay are subset to contain only the variable genes; default is TRUE}

\item{dims}{Which dimensions to use as input features, used only if features is NULL; defalut is 1:10}

\item{verbose}{Controls verbosity; default is TRUE}

\item{...}{Arguments passed to other methods}
}
\description{
Running analysis pipeline
}
