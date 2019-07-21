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
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A ggplot of cell and object projections showing the relative localization of the intracellular objects, grouped by cell length.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Renske van Raaphorst
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{createPlotList}}
}
\examples{

#get objects & meshes from WT D39 cells stained with Van-FL, analyzed with oufti.

library(bactMAP)

VanFL_mesh <- bactMAP::VanFL_meshWT
VanFL_obj <- bactMAP::VanFL_objWT

#plot objects

plotObjects(VanFL_obj, VanFL_mesh, groups=4, transparency=0.1)

#change colors/layout

library(ggdark)
library(ggplot2)
plotObjects(VanFL_obj, VanFL_mesh, groups=4, transparency=0.1, cellcolor="grey40", objectcolor="magenta") + ggdark::dark_theme_minimal()

}
