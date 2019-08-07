\name{extr_Original}
\alias{extr_OriginalStack}
\alias{extr_OriginalCells}
\alias{extr_Original}

\title{Convert original TIFF image stack used to gather spot and cell outline data.
}
\description{\code{extr_OriginalStack()} converts the original TIFF stack used to gather spot and cell outline data to a data frame containing the pixel x,y coordinates and the pixel value.

\code{extr_OriginalCells()} is connecting cell outlines to pixel values, using the output of \code{extr_OriginalStack()} and a \code{mesh} dataframe.
}
\usage{
extr_OriginalStack(picloc)

extr_OriginalCells(imdatframe, mesh, surroundings=FALSE, turnCell=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
extr_OriginalStack:
  \item{picloc}{file path of the .TIFF image stack.
}
extr_OriginalCells:
  \item{imdatframe}{output dataframe of extr_OriginalStack()}
  \item{mesh}{dataframe containing (at least) mesh coordinates, cell numbers and frame numbers.}
  \item{surroundings}{when TRUE, the function takes a box around the cells, when FALSE, only the pixels inside the cells are taken.}
  \item{turnCell}{when TRUE, the cells all turned in the same orientation. when FALSE, this is not the case.}
}
\details{
extr_OriginalStack needs the package "raster" to function.

extr_OriginalCells needs "SMDTools" and "raster" to function
}
\value{

The output of extr_OriginalStack is a list of dataframes including the x coordinates, y coordinates and intensity values of each pixel.

The output of extr_OriginalCells is a list of dataframes including:

\item{rawdata_turned}{a dataframe including the cell information and pixel coordinates per cell, all turned in the same direction.}
\item{mesh}{mesh data turned in the same direction}
}

\author{
Renske van Raaphorst
}
\note{

1. extr_OriginalCells is a heavy function since it's turning all pixels of each cell!

2. extr_OriginalStack works for single TIFFs and stacks, not hyperstacks. It only works with greyscale stacks.

}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
\dontrun{
#get tif
TIFFfilepath <- file.choose()

TIFFstack <- extr_OriginalStack(TIFFfilepath)

#plot stack
bactMAP::plotRaw(TIFFstack, frameN=1)

#get cell meshes
Meshpath <- file.choose()
meshes <- extr_Oufti(Meshpath)

#combine cells & tiff
extr_OriginalCells(TIFFstack, meshes$mesh)

}
}
