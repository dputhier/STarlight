% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/interactive.R
\name{fov_selection}
\alias{fov_selection}
\title{Interactive Visual Feature Picker}
\usage{
fov_selection(object = NULL)
}
\arguments{
\item{object}{An STGrid object. Alternatively, a dataframe or a matrix that contains
the coordinates (x, y) and the feature label. It thus should have at least three
columns x (x-coordinates), y (y-coordinates), and feature (feature labels).
Additional columns are supported.}
}
\value{
Launches a Shiny application which does not return a value to the R environment.
Selected features and configurations can be exported directly from the application interface.
}
\description{
This function launches an interactive Shiny application for visual exploration
and selection of features within a given dataset. It allows users to interactively
select features based on various criteria and visualization tools such as tangents,
areas, nearest points, and polygons.
}
\examples{
\dontrun{
example_dataset()
xen <- Xenium_Mouse_Brain_Coronal_7g
fov_selection(coord(xen))
data_set <- data.frame(x = runif(1000, 0, 1000), y = runif(1000, 0, 1000), feature = sample(c("Feature1", "Feature2"), 100, replace = TRUE))
fov_selection(data_set)
}
}
