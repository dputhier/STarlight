% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/st_grid_class.R
\name{as_matrix,STGrid-method}
\alias{as_matrix,STGrid-method}
\title{Convert the bin_mat slot of an STGrid object into an x/y matrix (e.g. to be displayed using image())}
\usage{
\S4method{as_matrix}{STGrid}(object, feature = NULL)
}
\arguments{
\item{object}{The STGrid object}

\item{feature}{The (meta)feature of interest.}
}
\description{
This function converts the bin_mat slot of an STGrid object into an x/y matrix (e.g. to be displayed using image())
}
\examples{
example_dataset()
xen <- Xenium_Mouse_Brain_Coronal_7g
image(as_matrix(xen))
}
