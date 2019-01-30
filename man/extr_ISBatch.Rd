\name{extr.ISBatch}
\alias{extr.ISBatch}
\alias{peakfitter}
\alias{ISBatch}

\title{
Upload ISBatch/Peakfitter data into R
}
\description{

}
\usage{
extr.ISBatch(dataloc)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataloc}{
  Path to the results table from the "peakfitter" option of the ISBatch imageJ plugin saved as comma-separated file.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Dataframe containing the following columns:
\item{$X1}{index number of the detected spot}
\item{$baseline}{baseline of peak intensity}
\item{$height}{height of detected peak}
\item{$x}{x coordinate of detected peak}
\item{$y}{y coordinate of detected peak}
\item{$sigma_x}{}
\item{$sigma_y}{}
\item{$fwhm_x}{full width half max (x direction) of detected peak}
\item{$fwhm_y}{full width half max (y direction) of detected peak}
\item{$fwhm}{average full width half max of detected peak}
\item{$error_baseline}{error of the baseline}
\item{$error_height}{error of the height}
\item{$error_x}{error of the x coordinate}
\item{$error_y}{error of the y coordinate}
\item{$error_sigma_x}{}
\item{$error_sigma_y}{}
\item{$error_fwhm_x}{}
\item{$error_fwhm_y}{}
\item{$error_fwhm}{}
\item{$frame}{frame of the image stack the spot was detected in}
}
\references{
\url{http://singlemolecule.github.io/iSBatch/}

Caldas, V.E., Punter, C.M., Ghodke, H., Robinson, A. and van Oijen, A.M., 2015. iSBatch: a batch-processing platform for data analysis and exploration of live-cell single-molecule microscopy images and other hierarchical datasets. \emph{Molecular BioSystems}, 11(10), pp.2699-2708.
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##Choose ISBatch peakfitter output file
dataloc <- file.choose()
##Convert to dataframe
spotframe <- extr.ISBatch(dataloc)
}

