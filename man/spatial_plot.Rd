% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/visualization.R
\name{spatial_plot}
\alias{spatial_plot}
\title{Plot x/y coordinates of molecules}
\usage{
spatial_plot(
  object = NULL,
  feat_list = NULL,
  colors = NULL,
  size = 0.1,
  coord_fixed = TRUE
)
}
\arguments{
\item{object}{The STGrid object.}

\item{feat_list}{The name of the features for which the spatial image will be created.}

\item{colors}{The colors to use for features in the spatial image.}

\item{size}{The size of the points}

\item{coord_fixed}{Logical value indicating whether to keep the aspect ratio fixed. Defaults to TRUE.}
}
\description{
Plot x/y coordinates of molecules of a spatial transcriptomics experiment.
}
\examples{
example_dataset()
xen <- Xenium_Mouse_Brain_Coronal_7g
spatial_plot(xen,
             feat_list=c("Chat", "Nwd2", "Ano1"),
             size=0.05)
}
\keyword{internal}
