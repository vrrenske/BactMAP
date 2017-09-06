\name{extr.ObjectJ}
\alias{extr.ObjectJ}
\alias{ObjectJ}

\title{Semi-manual entry of ObjectJ output data

}
\description{ObjectJ/ColiInspector/Chaintracer data input function. Semi-manual input.
}
\usage{
extr.ObjectJ(dataloc, spots, constr, xcolumns, ycolumns, xconstr, yconstr, spotcol,spotloc)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataloc}{path to tab delimeted .txt file containing cell dimension data and optionally manually-entered columns with spot x/y data and constriction x/y data}
  \item{spots="y"/"n"}{if \code{spots="y"}, the function will search for spot x/y data in the dataloc file. if spots is not specified and no seperate spot file is entered in the function, the function will prompt a question whether the user wants to include spots or not}
  \item{constr="y"/"n"}{as with \code{spots}, the function will ask the user whether the constriction x/y data is included in the \code{dataloc} file or not if not specified}
  \item{xcolumns}{If specified, name of manually defined column(s) containing x coordinates of fluorescent spots in the dataloc file. If there are more columns with the x coordinates of the spots, define the variable as c("column1", "column2", ...). If spots="n", the function will not ask for this variable. If the variable is not specified and spots="y", the function will ask for this variable while running.}
  \item{ycolumns}{If specified, name of manually defined column containing y coordinates of fluorescent spots in the dataloc file. If there are more columns with y coordinates of the spots, define the variable as c("column1", "column2", ..). If spots="n", the function will not ask for this variable. If the variable is not specified and spots="y", the function will ask for this variable while running.}
  \item{xconstr}{If specified, name of manually defined column containing x coordinates of the constriction points in the dataloc file. If constr="n", the function will not ask for this variable. If the variable is not specified and constr="y", the function will ask for this variable while running.}
  \item{yconstr}{If specified, name of manually defined column containing y coordinates of the constriction points in the dataloc file. If constr="n", the function will not ask for this variable. If the variable is not specified and constr="y", the function will ask for this variable while running.}
  \item{spotcol}{If spots are detected, it is possible to add the name of the column containing the number of spots/cell. If not defined, the function will ask for it while running.}
  \item{spotloc}{path to the comma-separated file generated when "calculate spot positions" is runned in ObjectJ/ColiInspector. When spotloc is specified, spots, xcolumns and ycolums should not be used.}
}

\details{
The only obliged input variable is "dataloc", containing the cell length and width calculated by ObjectJ/ColiInspector. It is optional to include spot localizations in the dataloc file and specify them either in the function call, or the function will ask for it while running. It is also possible to add a ObjectJ file with the relative location of the spots inside the cell.
To get the ObjectJ file with the relative location of the spots inside the cell, use the function "calculate spot positions" in ObjectJ to generate a seperate spot-location file. Otherwise, it is also possible to calculate the spot positions detected in ObjectJ or another spot detection program relative to the cell dimensions using the function \code{shinyspots::putSpotsinBox()}
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Dataframe containing the length and width of the cells, and if selected, spot and/or constriction point coordinates.
\item{$cell}{number identifying the bacterial cell}
\item{$length}{length (pixel) of the longest axis of the cell}
\item{$max.width}{the diameter, or maximum width of the cell in pixels}
\item{$frame}{the nth frame of the image stack}

When spot x/y coordinates are indicated in the dataloc file:
\item{$x}{x coordinate of the indicated spot}
\item{$y}{y coordinate of the indicated spot}
\item{$spotnum}{number of spots per cell}

When constriction point distances to the length/width axis are indicated in the dataloc file:
\item{$Lc}{distance of constriction point to the length axis}
\item{$Ld}{distance of constriction point to the width axis}


When calculated spots are added as a seperate results table:
\item{$L}{The shortest distance of the spot to the length axis}
\item{$D}{The shortest distance of the spot to the width axis}
}

\references{
Vischer, N.O., Verheul, J., Postma, M., van Saparoea, B.V.D.B., Galli, E., Natale, P., Gerdes, K., Luirink, J., Vollmer, W., Vicente, M. and den Blaauwen, T., 2015. Cell age dependent concentration of \emph{Escherichia coli} divisome proteins analyzed with ImageJ and ObjectJ. \emph{Frontiers in microbiology}, 6.

Syvertsson, S., Vischer, N.O., Gao, Y. and Hamoen, L.W., 2016. When Phase Contrast Fails: ChainTracer and NucTracer, Two ImageJ Methods for Semi-Automated Single Cell Analysis Using Membrane or DNA Staining. \emph{PloS one}, 11(3), p.e0151267.
}
\author{
Renske van Raaphorst
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x)
{
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
