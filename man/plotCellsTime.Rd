\name{plotCellsTime}
\alias{plotCellsTime}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Function to plot cell fluorescence from a TIFF over time.
}
\description{
  This function creates a cell tower or movie based on the pixel values and cell outlines of the cell. As input, it takes the 'rawdata_turned' dataframe output from \code{\link{extr_OriginalCells}}.
}
\usage{
plotCellsTime(celdat,
              updown = T,
              movie = F,
              viridisoption = "magma",
              cellN,
              minf,
              maxf,
              outlines=FALSE,
              meshdata)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{celdat}{
  TIFF and mesh outline data from a timelapse movie as it was obtained by running \code{extr_OriginalCells}.
}
  \item{updown}{
  When \code{updown}==\code{TRUE}, a horizontal cell tower is made, otherwise, the cells are displayed left-right. \code{movie}=={TRUE} will overwrite this and return an animation instead.
}
  \item{movie}{
When \code{movie} == \code{TRUE}, an animation made with \href{https://github.com/thomasp85/gganimate}{gganimate} is returned.
}
  \item{viridisoption}{
  The color palette used. The default palettes are part of the \href{https://cran.r-project.org/web/packages/viridis/vignettes/intro-to-viridis.html}{viridis} color palettes (default = "magma").
}
  \item{cellN}{
  The cell number or numbers to be displayed. When not indicated, a list of plots of all cells will be made.
}
  \item{minf}{
  The minimum frame to be displayed.
}
  \item{maxf}{
  The maximum frame to be displayed
}
  \item{outlines}{when TRUE and \code{meshdata} is specified, the cell outlines are drawn over the pixels}
  \item{meshdata}{to be able to draw outlines, specify your mesh dataframe.}
}

\value{
A plot, list of plots or animation of a single cell's fluorescence over time.
}

\author{
Renske van Raaphorst
}

\examples{
\dontrun{
##In the BactMAP example data, the raw data of one cell is saved as "TurnedCell4":
##Download this via https://veeninglab.com/f/example_datasets.zip

load("TurnedCell4.rda")
#to plot a left-right tower:
plotCellsTime(TurnedCell4, updown=F, viridisoption="inferno")

#or to create an animation:
plotCellsTime(TurnedCell4, movie=T)
}
}
