\name{extr.Original}
\alias{extr.OriginalStack}
\alias{extr.OriginalCells}

\title{Convert original TIFF image stack used to gather spot and cell outline data.
}
\description{extr.OriginalStack() converts the original TIFF stack used to gather spot and cell outline data to a data frame containing the pixel x,y coordinates and the pixel value.

extr.OriginalCells() is connecting the cells from a
}
\usage{
extr.OriginalStack(picloc)

extr.OriginalCells(imdatframe, mesh)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
extr.OriginalStack:
  \item{picloc}{file path of the .TIFF image stack.
}
extr.OriginalCells:
  \item{imdatframe}{output dataframe of extr.OriginalStack()}
  \item{mesh}{dataframe containing (at least) mesh coordinates, cell numbers and frame numbers.}
}
\details{
extr.OriginalStack needs the package "raster" to function.

extr.OriginalCells needs "SMDTools" and "raster" to function
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
extr.OriginalCells is a heavy function since it's turning all pixels of each cell!
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
