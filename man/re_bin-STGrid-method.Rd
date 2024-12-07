% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/st_grid_class.R
\name{re_bin,STGrid-method}
\alias{re_bin,STGrid-method}
\title{Re-bin a STGrid object.}
\usage{
\S4method{re_bin}{STGrid}(object, bin_size, verbose = TRUE)
}
\arguments{
\item{object}{The STGrid object.}

\item{bin_size}{The size of the bin.}

\item{verbose}{Whether to be display progress bar.}
}
\description{
Re-bin a STGrid object.
}
\examples{
example_dataset()
res <- re_bin(Xenium_Mouse_Brain_Coronal_7g, bin_size=10)
}
