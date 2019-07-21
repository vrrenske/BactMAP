\name{bactKymo}
\alias{bactKymo}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Creates kymographs or demographs from pixel values & cell dimensions.
}
\description{

bactKymo() takes the output of \code{extr_OriginalCells()} and makes kymographs or demographs of them, depending on the chosen options. By default, all cells are used to make a demograph where the average fluorescence over the length axis of each cell is plotted, where all cells are ordered by cell length.

When 'timeD' is set to TRUE, bactKymo() makes kymographs of single cells. However, when 'percDiv' is set to true as well, bactKymo will make a kymograph of the average of all cell's fluorescence over their length axis, binning cells by the moment of division they are in (see \code{perc_Division()}).

Other settings which can be changes are the bin size (in how many groups is the average fluorescence intensity over the cell length calculated), whether the bins are placed such, that they represent the cell length at that moment (sizeAV=TRUE, takes more computing time), and when using a demograph, whether the extreme values are removed (cutoff_demograph).

The output will standardly be in the "viridis" color scale of the viridis package. However, this can be overwritten manually by adding a fill color scale to the plot.

To make computing time shorter while playing with the settings, it is possible to first prepare the kymographs using \code{prepForKymo()}. This is generally not necessary, but recommended when plotting single cell kymographs one by one.

}
\usage{
bactKymo(originalCells,
          timeD = FALSE,
          dimension = "length",
          bins = 25,
          sizeAV = FALSE,
          cells = "all",
          prep = TRUE,
          percDiv = FALSE,
          cutoff_demograph = 0.975,
          mag,
          legend=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{originalCells}{
  The output of extr_OriginalCells, which combines TIFF data with mesh data.
}
  \item{timeD}{
  if \code{timeD==TRUE}, the resulting kymographs will be sorted by time, not by cell length and will be made of single cells.
}
  \item{dimension}{
  default="length", but when put to "width" the kymographs will average over cell width.
}
  \item{bins}{
  the amount of bins the fluorescence intensity will be averaged over (per cell/time point)
}
  \item{sizeAV}{
  when \code{sizeAV==TRUE}, the kymographs will reflect the cell length (or width) on the y axes. Default to FALSE because the computing time is faster when there is no sizing.
}
  \item{cells}{
  default = "all". 'cells' is a numeric vector OR the word "all". When set  to "all" ,all cells will be processed into a demograph (when timeD==FALSE) or a kymograph will be made of all cells (when timeD==TRUE). When a numeric vector of length 1 or more, the cells identified with the cell numbers in these vector will be processed into a demograph (when timeD == FALSE), or a kymograph will be made of all single cells (when timeD==TRUE). Finally, when only 1 cell number is identified, a kymograph will be made of this cell (when timeD == TRUE). It is not possible to make a plot of a single cell without time dimension.
}
  \item{prep}{
  default = TRUE. Inside the function, the data will be prepared so it is possible to plot kymographs. Sometimes, when you want to run \code{bactKymo()} multiple times on the same dataset, it can save time to prepare the dataset only once. For this you can use the function \code{prepForKymo()} first and use the as the argument \code{originalCells}, while setting 'prep' to FALSE.
}
  \item{percDiv}{
  this will call for \code{perc_Division()} and bin the average fluorescence over cell length (or width) in 10 groups, from just divided to almost dividing. It will, by default, plot all cells and only works with timelapse data.
}
  \item{cutoff_demograph}{
  to make it easier to look at demographs of all cells, the highest 2.5 percent of the intensity values are removed from the demograph. set the cutoff to 1 to not remove anything, or put it lower to remove more values.
}
  \item{mag}{
  magnification converter (see \code{\link{Pixels2um}}. If used, the pixel values will be converted to micron.
}
  \item{legend}{
  when \code{legend==TRUE}, a legend will be shown on the right. when FALSE, the legend will be omitted.
  }
}

\details{
Many options are possible within \code{bactKymo()}. Make sure you set the time dimension & cells correctly! For examples, see \href{https://github.com/vrrenske/BactMAP/wiki}{the BactMAP wiki}

Note that when you have an output of many plots, it is useful to use \code{gridExtra::grid.arrange()} to save all plots in a PDF, this makes the plots easier to investigate and saves loading time in between.

}
\value{
One or multiple demographs/kymographs (see the details in \code{arguments}).
}

\author{
Renske van Raaphorst
}

\examples{

#to plot a timelapse demograph of 1 cell:

kymo_cell4 <- bactKymo(originalCells = bactMAP::TurnedCell4, timeD=TRUE, sizeAV=TRUE)

#to change the color scheme:
OrangeHot <- getPalette("OrangeHot")
kymo_cell4 + scale_fill_gradient2(low=OrangeHot[1], mid=OrangeHot[2], high=OrangeHot[3])


}
