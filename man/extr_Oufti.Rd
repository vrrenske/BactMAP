\name{extr.Oufti}
\alias{extr.Oufti}
\alias{Oufti}
\alias{extr.OuftiCSV}
\alias{OuftiCSV}

\title{Function to upload Oufti meshes and spot localizations directly from .Mat files to R.
}
\description{

}
\usage{
extr.Oufti(matfile)

extr.OuftiCSV(csvfile)
}

\arguments{
  \item{matfile}{.MAT output of oufti segmentation. whether spot data is included will automatically be detected. not possible to include object data for now.}
  \item{csvfile}{.CSV or .TXT output of oufti segmentation (saved with "tools"-->"export to csv" in the Oufti GUI). Whether spot data and/or object data is included will automatically be detected.}
%%     ~~Describe \code{x} here~~
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A list containing:
\item{$cellList}{The complete cellList as it is displayed in oufti's matlab output.}
\item{$mesh}{A dataframe containing the following columns: \code{x, y, cell, frame, num, max.width, length, steplength, max_length}. See standard mesh output for more information.}
\item{$spots}{A dataframe containing the following columns: \code{l, d, x, y, position, adj_Rsquared, CI_xy, frame, cell}. See \url{www.oufti.org} for more information on the spot data output.}

optional in case of \code{extr.Oufti}:
\item{$obj}{A dataframe containing the coordinates of detected objects.}
}

\references{
\url{www.oufti.org}

Paintdakhi A, et al. (2016) Oufti: an integrated software package for high-accuracy, high-throughput quantitative microscopy analysis. \emph{Molecular Microbiology} 99(4):767-777.
}
\author{
Renske van Raaphorst
}
\note{

}



\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##upload matlab oufti file
matfile <- file.choose()
##extract data using extr.Oufti()
output <- extr.Oufti(matfile)
##get spot dataframe
spots <- output$spots
##plot spots of frame no.1
ggplot2::ggplot(spots[spots$frame==1,], ggplot2::aes(x=x, y=y, color=as.factor(cell))) + ggplot2::geom_point() + ggplot2::theme_minimal()
}

