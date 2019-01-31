\name{extr_SuperSeggerCells}
\alias{extr_SuperSeggerCells}
\alias{extr_SuperSegger}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Convert SuperSegger cell outlines to standard BactMAP format
}
\description{
Part of the *extr_*-functions, extr_SuperSeggerCells takes the mask created by SuperSegger of each cell (and if available, timepoint) and saves this as a *mesh* dataframe.
}
\usage{
extr_SuperSeggerCells(loc, frames, mag)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{loc}{
The file path where SuperSegger's output files are located
}
  \item{frames}{
The number of xy frames (so not timepoints, just different locations!) imaged
}
  \item{mag}{
magnification conversion factor name (which is part of Pixels2um)
}
}
\details{
Use *addPixels2um()* to add a new conversion factor.
}
\value{
 \item{cellList}{dataframe with content similar to SuperSegger's "Clist"}
 \item{mesh}{dataframe containing cell coordinates and dimensions}
 \item{pixel2um}{the magnification conversion factor used}
}
\references{
Stylianidou, Stella, et al. "SuperSegger: robust image segmentation, analysis and lineage tracking of bacterial cells." Molecular microbiology 102.4 (2016): 690-700.
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


