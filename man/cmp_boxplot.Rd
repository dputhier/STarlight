% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/visualization.R
\name{cmp_boxplot}
\alias{cmp_boxplot}
\title{Create a boxplot/jitter plot to show molecule counts distribution.}
\usage{
cmp_boxplot(
  object,
  normalized = TRUE,
  transform = c("None", "log2", "log10", "log"),
  colors = c("#3074BB", "#BE5B52"),
  ...
)
}
\arguments{
\item{object}{A STCompR object.}

\item{normalized}{A logical value indicating whether to use normalized counts.}

\item{transform}{A character string specifying the transformation to be applied to molecule counts.
Options include "None", "log2", "log10", and "log".}

\item{colors}{A vector of colors for the boxplot/jitter plot.}

\item{...}{Additional arguments to be passed to the underlying ggpol::geom_boxjitter function.}
}
\value{
A ggplot object displaying the molecule counts distribution.
}
\description{
Create a boxplot/jitter plot to show molecule counts distribution.
}
\examples{
example_dataset("11284296/files/cmp_xen")
cmp_xen
cmp_boxplot(cmp_xen, normalized = TRUE, transform = "log2", colors = c("blue", "red"))

}
\keyword{internal}
