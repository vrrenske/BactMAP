\name{prepForKymo}
\alias{prepForKymo}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Prepare whole dataset for single cell kymographs.
}
\description{
Extension from \code{bactKymo}, this function prepares a dataset for kymograph plotting. Recommended if you plan to plot a lot of kymographs from one dataset one by one;  this will speed up the process. If you don't plan to do this, use \code{bactKymo} immediately.
}
\usage{
prepForKymo(turnedCells, dimension = "length", bins = 25, sizeAV = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{turnedCells}{
  Output from \code{extr_OriginalCells}, which includes the pixel values (of the fluorescence signal) in each cell.
}
  \item{dimension}{
 Default = "length", "width" is also possible. This decides on which dimension the kymograph will be calculated.
}
  \item{bins}{
 Default = 25. The amount of bins in which the average fluorescence over length or width will be calculated.
}
  \item{sizeAV}{
 When \code{sizeAV == TRUE}, the bins are placed along the y axis of the plot according to cell size. This takes more computation time.
}
}

\value{
A dataframe suitable for \code{bactKymo()}.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{



}
