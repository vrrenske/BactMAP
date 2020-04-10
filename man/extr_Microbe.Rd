\name{extr_MicrobeJ}
\alias{extr_MicrobeJ}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{upload MicrobeJ output into R for further analysis
%%  ~~function to do ... ~~
}
\description{
Upload the output of MicrobeJ into R. Upto now the code supports fluorescent localization data & cell outlines, both saved as seperate .CSV or .TXT files in the MicrobeJ GUI by first pressing \emph{<right click> -> expand all} on one of the Results columns and then pressing \emph{<right click> -> export} for both the cell contour dataset and the fluorescent localization dataset.
}
\usage{
extr_MicrobeJ(dataloc, spotloc, objectloc, mag = "No_PixelCorrection", sepspot=",", sepmesh=",", sepobj=",", cellList=FALSE, keeprealvalues=FALSE, magcor = c("dataloc", "spotloc", "objectloc"))
}


\arguments{
  \item{dataloc}{path of the exported microbeJ file (.csv or .txt) containing the outlines of the bacterial cells.}
  \item{spotloc}{path of the exported microbeJ file (.csv or .txt) containing the sub-directory with fluorescent spot coordinates}
  \item{objectloc}{path of the exported microbeJ file (.csv or .txt) containing the outlines of objects/filaments in the cell.}
  \item{mag}{the pixel conversion factor to be used, default = "No_PixelCorrection"}
  \item{sepspot}{default is ",". when the to-be-imported .csv/txt is seperated by another value, change here.}
  \item{sepmesh}{default is ",". when the to-be-imported .csv/txt is seperated by another value, change here.}
  \item{sepobj}{default is ",". when the to-be-imported .csv/txt is seperated by another value, change here.}
  \item{cellList}{when cellList==TRUE, the dataframe 'cellList' (see Values below) will be in the output list.}
  \item{keeprealvalues}{default FALSE. BactMAP checks whether the mesh and spot files are both in pixels or in micron and corrects this. When keeprealvalues==TRUE, BactMAP skips this step.}
  \item{magcor}{indicates the input datasets which need to be corrected (pixel to micron). if not changed, all datasets are corrected.}
  \item{extracols}{if there are any other main variables inside 'dataloc' you wish to pass on to your output, add them here as a vector of character strings (for instance "myvariable" for 1 variable, or c("myvariable1", "myvariable2") for 2 variables.)}

}

\details{
Input either \code{dataloc}, \code{spotloc} or both to receive a list of dataframes containing the microbeJ data.

Use \code{getPixels2um()} and \code{addPixels2um()} to check which conversion factors are loaded and to add new ones.

}
\value{
List of items containing:
\item{$mesh}{dataframe containing MicrobeJ outlines of the cell}
\item{$spotframe}{dataframe containing MicrobeJ spot x/y coordinates and cell identifying number.}
\item{$spots_relative}{dataframe containing the relative x/y coordinates of the spots inside the cell}
\item{$cellList}{(a list of) the original input dataframe(s)}
\item{$pixels2um}{the pixel conversion factor used in the function}
}


\references{
\url{https://microbej.com}

Ducret, A., Quardokus, E.M. and Brun, Y.V., 2016. MicrobeJ, a tool for high throughput bacterial cell detection and quantitative analysis. \emph{Nature microbiology}, 1, p.16077.

}
\author{
Renske van Raaphorst
}

\examples{
\dontrun{
##upload MicrobeJ datasets
dataloc <- file.choose()
spotloc <- file.choose()
##extract data
micdata <- extr_MicrobeJ(dataloc, spotloc, mag = "100x_LeicaVeening")

}
}

