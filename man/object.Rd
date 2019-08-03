\name{object}
\alias{object}
\docType{data}
\title{
Object DataFrames
}
\description{
A dataframe containing object information, either \code{objectframe} or \code{object_relative}.

These dataframes can be used in various of the bactMAP functions: \code{\link{plotTracks}, \link{plotObjects}, \link{plotOverlay} & \link{plotRaw}}.
}
\format{
  In the case of \code{objectframe}, the most basic of the two, the object coordinates and cell information are given:
  \describe{
    \item{\code{ob_x}}{the x-coordinate (in pixels) of the object shape outline}
    \item{\code{ob_y}}{the y-coordinate (in pixels) of the object shape outline}
    \item{\code{obnum}}{the number of the object inside one cell}
    \item{\code{obpath}}{number indicates the order of the object coordinates, so a polygon will be drawn correctly}
    \item{\code{frame}}{the number of the image frame}
    \item{\code{cell}}{the number of the bacterial cell}
    \item{\code{obID}}{unique object identifier}
  }

  In the case of \code{object_relative}, the following colums are added:
  \describe{
  \item{\code{max.length}}{the cell length of the given cell the object is in.}
  \item{\code{num}}{a rank for the cell size - the smaller \code{num}, the smaller the cell size}
  \item{\code{max.width}}{the cell width of the given cell this object is in.}
  \item{\code{Dum}}{the relative localization of the mid-point of the object on the width axis of the cell, in micron.}
  \item{\code{Lmid}}{the relative localization of the mid-point of the object on the length axis of the cell, in micron.}
  \item{\code{ob_out_x}}{the relative x-coordinate (in micron) of the object shape outline to the cell}
  \item{\code{ob_out_y}}{the relative y-coordinate (in micron) of the object shape outline to the cell}
  \item{\code{pole1}}{the distance of one cell pole to mid-cell}
  \item{\code{pole2}}{the distance of the other cell pole to mid-cell. Since mid-cell is placed at [0,0], this is always the negative of pole1.}
  \item{\code{max_um}}{\code{max.length} in micron}
  \item{\code{maxwum}}{\code{max.width} in micron}

  }
}
\details{
An example of an object dataframe can be taken from the bactMAP package. How is shown in "examples" below.
}

\examples{
\dontrun{
#take VanFL object data from https://veeninglab.com/f/example_datasets.zip
load("VanFL_objWT.rda")

summary(VanFL_objWT.rda)

}
}

\keyword{datasets}
