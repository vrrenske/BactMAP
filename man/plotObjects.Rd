\name{plotObjects}
\alias{plotObjects}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Create plots similar to \link{createPlotList} using intracellular object outlines instead of points.
}
\description{
Plot projection of fluorescent objects inside cells in one or more cell groups based on cell length.
}
\usage{
plotObjects(obdat,
            meshdat,
            groups = 1,
            cellcolor = "black",
            objectcolor = "green",
            transparency = 0.1,
            getdata = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obdat}{
  Object dataset "object_relative" as it is extracted using one of the "extr_"-functions.
}
  \item{meshdat}{
  A mesh dataset as it is extracted using one of the "extr_"-functions (optional).
}
  \item{groups}{
  The amount of groups to devide the cell projections in (by cell length). Default = 1.
}
  \item{cellcolor}{
  The color of the cells in the cell projections. Default = black.
}
  \item{objectcolor}{
  The color of the objects in the projections. Default = green.
}
  \item{transparency}{
  The transparency of the objects and cells in the overlay. Default = 0.1. When analyzing ~1000 cells, it is advised to set this number to ~0.02 or less.
}
  \item{getdata}{
  When getdata==TRUE, the function returns a dataset with the grouped cells next to the plot. Default is set to FALSE.
}

}

\value{
A ggplot of cell and object projections showing the relative localization of the intracellular objects, grouped by cell length.
}

\author{
Renske van Raaphorst
}

\seealso{
\code{\link{createPlotList}}
}
\examples{
\dontrun{
#get objects & meshes from WT D39 cells stained with Van-FL, analyzed with oufti.
#download from https://veeninglab.com/f/example_datasets.zip

load("VanFL_meshWT")
load("VanFL_objWT")

#plot objects

plotObjects(VanFL_objWT, VanFL_meshWT, groups=4, transparency=0.1)
}
}
