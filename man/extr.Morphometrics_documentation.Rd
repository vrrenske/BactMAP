\name{extr.Morphometrics}
\alias{extr.Morphometrics}
\alias{Morphometrics}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Import Morphometrics .MAT CellList Output Into R
%%  ~~function to do ... ~~
}
\description{Function returns the CellList out of the Oufti .MAT output as a dataframe.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
extr.Morphometrics(morphpath)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{morphpath}{path of the .MAT file containing the Oufti output
%%     ~~Describe \code{x} here~~
}
}
\details{
This function works with output from morphometrics, a bacterial cell segmentation program for phase contrast images which works in Matlab. It filters the information in the output of morphometrics to the values below. If you are interested in the other values in the .MAT output, use \emph{extr.Morphometrics.cellList} to upload the full file.
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Dataframe containing:
\item{Xcont}{The X coordinates of the cell contour point}
\item{Ycont}{The Y coordinates of the cell contour point}
\item{CellID}{The number identifying the bacterial cell}
\item{Area}{The area of the bacteria cell in pixel^2}
\item{poleL}{The point "numpoint" where "Xcont" and "Ycont" represent the x,y coordinates of the left pole of the bacterial cell}
\item{poleR}{The point "numpoint" where "Xcont" and "Ycont" represent the x,y coordinates of the right pole of the bacterial cell}
\item{frame}{The number identifying the nth frame of the image stack}
\item{numpoint}{The nth point of the coordinates of the outline of the bacterial cell}

%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
\url{https://simtk.org/projects/morphometrics}

Ursell, T., Lee, T.K., Shiomi, D., Shi, H., Tropini, C., Monds, R.D., Colavin, A., Billings, G., Bhaya-Grossman, I., Broxton, M. and Huang, B.E., 2017. Rapid, precise quantification of bacterial cellular dimensions across a genomic-scale knockout library. \emph{BMC biology}, 15(1), p.17.
}
\author{
Renske van Raaphorst
%%  ~~who you are~~
}
\note{
This code uses the R.matlab package: \url{https://cran.r-project.org/web/packages/R.matlab/R.matlab.pdf}
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{extr.Morphometrics.cellList}}
}
\examples{
##Choose Morphometrics output file
morphpath <- file.choose()
##Get dataframe containing cell outlines
morph.outlines <- shinyspots::extr.Morphometrics(morphpath)
##Plot all cell outlines
ggplot2::ggplot(morph.outlines, ggplot2::aes(x=Xcont, y=Ycont)) + ggplot2::geom_point()
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{IO}% use one of  RShowDoc("KEYWORDS")
\keyword{file}% __ONLY ONE__ keyword per line
