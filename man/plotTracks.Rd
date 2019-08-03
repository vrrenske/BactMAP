\name{plotTracks}
\alias{plotTracks}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot tracked spots and/or objects over time in 1 or 2 dimensions.
}
\description{
With this function, you can plot spot and/or objects over time, either as a plot over time (x axis) versus localiton on the length or width axis of the cell, or as a cell projection, where color indicates time. When only one spatial dimension is displayed, it is possible to plot trajectories of different fluorescence channels. When the function detects that different conditions are used, the plots will be faceted by condition.
}
\usage{
plotTracks( meshdata,
            spotdata,
            objectdata,
            tracks = FALSE,
            ignore_singles = FALSE,
            movie = FALSE,
            timepalette_lines = "viridis",
            timepalette_fill = "magma",
            dimension = c("length", "width"),
            turn_cells = TRUE,
            mag,
            cell = "all",
            transparency = 0.2
            )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{meshdata}{
  A \code{\link{mesh}} dataframe containing cell outlines.
}
  \item{spotdata}{
  A \code{\link{Spot Data}} dataframe containing spot localizations.
}
  \item{objectdata}{
  A \code{\link{object}} dataframe containing fluorescent object outlines.
}
  \item{tracks}{
  When \code{tracks} == \code{TRUE}, the track of individual spots over time will be shown. When \code{tracks} == \code{FALSE}, only the spots will be shown.
}
  \item{ignore_singles}{
  When \code{ignore_singles} == \code{FALSE}, all detected spots will be shown. When \code{ignore_singles} ==\code{TRUE}, only those spots will be shown which are tracked over time.
}
  \item{movie}{
  When \code{movie}==\code{FALSE} (default), the output will be a static plot. When \code{movie} == \code{TRUE},  \href{https://github.com/thomasp85/gganimate}{gganimate} (version => 1.0.3) will be used to make a GIF animation.
}
  \item{timepalette_lines}{
  A \href{https://cran.r-project.org/web/packages/viridis/vignettes/intro-to-viridis.html}{viridis} palette for both the spots & tracks as the outlines of the meshes, where color indicates the time in frames. Default = \code{"viridis"}.
}
  \item{timepalette_fill}{
A \href{https://cran.r-project.org/web/packages/viridis/vignettes/intro-to-viridis.html}{viridis} palette for the fill of the objects, where color indicates the time in frames. Default = \code{"magma"}.
}
  \item{dimension}{
Either  \code{"length"},  \code{"width"}, or  \code{c("length", "width")}(default). Indicates the dimensions which will be plotted next to time (frames).
}
  \item{turn_cells}{
When cell projections are plotted (both cell length & width), it is possible to plot cells with their original x/y coordinates, or turn them so their length axis is parallel to the x-axis with mid-cell at 0. When plotting only one cell dimension, this command is ignored and cells are turned by default.
}
  \item{mag}{
The magnification converter (see \code{\link{pixels2um}}).
}
  \item{cell}{
The cell number or numbers as a numeric vector (e.g. \code{c(1,2,10)}), or (default) "all". When more then one cell are indicated, the function returns a plot faceted by cell number.
}
  \item{transparency}{
  Transparency of both the objects as the cell polygons. The cell polygons are set to be 5 times more transparent than the object polygons.
}
}

\value{
A plot or a faceted plot of the trajectory of the cell shape, fluorescent spots and/or the localization of detected fluorescent objects in the cell over time. The plots can either by location on the length-axis over time, location on the width-axis over time, or x/y location and shape over time.
}

\author{
Renske van Raaphorst
}

\examples{
\dontrun{
  ##Get example data from bactMAP on https://veeninglab.com/f/example_datasets.zip

  ##tracked FtsZ
  load("FtsZ_tracks.rda")

  ##and corresponding meshes
  load("DnaX_mesh.rda")

  #plot ftsZ on length axis/time
  plotTracks(meshdata=meshes,
              spotdata=FtsZ_tracks,
              dimension="length",
              mag="100x_DVMolgen",
              tracks=TRUE,
              cell=4)
}

}
