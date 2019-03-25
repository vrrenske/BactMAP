\name{extr_Original}
\alias{extr_Original}
\alias{extr_OriginalStack}
\alias{extr_OriginalCells}

\title{Convert original TIFF image stack used to gather spot and cell outline data.
}
\description{extr_OriginalStack() converts the original TIFF stack used to gather spot and cell outline data to a data frame containing the pixel x,y coordinates and the pixel value.

extr.OriginalCells() is connecting the cells from a
}
\usage{
extr_OriginalStack(picloc)

extr_OriginalCells(imdatframe, mesh)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
extr_OriginalStack:
  \item{picloc}{file path of the .TIFF image stack.
}
extr_OriginalCells:
  \item{imdatframe}{output dataframe of extr_OriginalStack()}
  \item{mesh}{dataframe containing (at least) mesh coordinates, cell numbers and frame numbers.}
}
\details{
extr_OriginalStack needs the package "raster" to function.

extr_OriginalCells needs "SMDTools" and "raster" to function
}
\value{

}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{

1. extr_OriginalCells is a heavy function since it's turning all pixels of each cell!

2. extr_OriginalStack works for single TIFFs and stacks. It only works with greyscale stacks.

}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

#get tif
TIFFfilepath <- system.file("extdata", "DnaX-GFP.TIF", package = "bactMAP"))

TIFFstack <- extr_OriginalStack(TIFFfilepath)

#plot stack
bactMAP::plotRaw(TIFFstack, frameN=1)

#get cell meshes
Meshpath <- system.file("extdata", "DnaX-mesh.mat", package="bactMAP")
meshes <- extr_Oufti(Meshpath)

#combine cells & tiff
extr_OriginalCells(TIFFstack, meshes$mesh)


}
