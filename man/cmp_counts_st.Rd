% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/visualization.R
\name{cmp_counts_st}
\alias{cmp_counts_st}
\title{Compare feature counts across multiple STGrid objects.}
\usage{
cmp_counts_st(
  ...,
  features = NULL,
  type = c("barplot"),
  names = NULL,
  transform = c("None", "log2", "log10", "log"),
  fill_color = NULL,
  border_color = "black"
)
}
\arguments{
\item{...}{A set of STGrid objects.}

\item{features}{A character vector specifying the features (e.g., genes) to compare. Default is NULL.}

\item{type}{Type of plot to generate. Currently only barplot is supported.}

\item{names}{Optional character vector specifying names for each STGrid object. If NULL, default names will be assigned.}

\item{transform}{Transformation method for the counts. Options are "None" (default), "log2", "log10", or "log".}

\item{fill_color}{Optional character vector specifying fill colors for different conditions in the plot.}

\item{border_color}{Color for the border of bars in bar plots. Default is "black".}
}
\value{
A ggplot object representing the comparison of counts across different conditions or samples.
}
\description{
This function compares feature counts across multiple STGrid objects.
}
\examples{
example_dataset()
xen <- Xenium_Mouse_Brain_Coronal_7g
x_bins <-  bin_x(xen)[181:nbin_x(xen)]
y_bins <-  bin_y(xen)[101:nbin_y(xen)]
xen_r1 <- xen[x_bins, y_bins]
x_bins <-  bin_x(xen)[61:101]
y_bins <-  bin_y(xen)[101:nbin_y(xen)]
xen_r2 <- xen[x_bins, y_bins]

cmp_counts_st(xen, xen_r1, xen_r2, features = c("Chat", "Ano1"), normalized = FALSE)
cmp_counts_st(xen, xen_r1, xen_r2,
             features = c("Chat", "Ano1"),
              normalized = TRUE, fill_color=grDevices::rainbow(3))
}
