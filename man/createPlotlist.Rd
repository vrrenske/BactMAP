\name{createPlotlist}
\alias{createPlotlist}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to create descriptive plots of spot localization per group of cells
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
createplotlist <- function(Spots, inp=4, bycol="length", MESH, colorpalette="GreenYellow", mag="100x_LeicaVeening", AllPlot=T, Xm="X", Ym="Y")
}
\arguments{
  \item{Spots}{Data frame containing -at least- column with relative spot x-localizations, y-localizations, cell ID and frame number.}
  \item{inp}{Number of groups the data will be divided in. Default = 4.}
  \item{bycol}{Identifies the column based on which the data will be divided. Default = "length".}
  \item{MESH}{Data frame containing -at least- the x and y coordinates of the cell meshes, the maximum cell width, maximum cell length, cell ID and frame.Optional.}
  \item{colorpalette}{Color palette the heatmap plots will be plotted in. Default is "GreenYellow". Check which palettes are loaded by \code{showColorPaletes()} and add your own palette using \code{addPallete()}}
  \item{mag}{The name of the conversion factor for pixel to micrometer. Default is "100x_LeicaVeening" which is 0.0499538. Check \code{\link{getPixels2um}} for the loaded conversionfactors and use \code{\link{addPixels2um("ConversionName", conversionfactor)}} to add a conversion factor to the list.}
  \item{AllPlot}{When \code{TRUE}, the output x/y plots will have sideplots with histograms of the spot localization. When \code{FALSE}, only the x/y plots will be plotted.}
  \item{Xm and Ym}{By default, the column names of the x/y coordinates in the MESH data frames are "X" and "Y". When this is not the case, one can tell the function by using \code{Xm} and \code{Ym}. Enter a character string with the name of the x-coordinate column for Xm, and the name of the y-coordinate column for Ym.}


}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A list consisting of the following items:
\item{$lengthplot}{Heatmap histogram plot of the length-axis localizations of the spots ordered by the length of the cell.}
\item{$widthplot}{Heatmap histogram plot of the width-axis localizations of the spots ordered by the length of the cell.}
\item{$qplots}{If the outlines of the cell are not added to the function (dataframe \code{MESH} not indicated), the output will be a list of plots of length \code{inp}, with relative x/y coordinates of the spots plotted as a density function. The data is grouped by the variable \code{bycol}, which is \code{length} by default.If the outlines of the cells are added to the function (dataframe \code{MESH} indicated), the average cell outline of each group will be overlayed on top of each plot. If \code{AllPlot==TRUE}, histograms of the x and y coordinates of the spots will be added as sideplots.}
\item{$plottotal}{Relative x/y coordinates of all spots with average cell outline overlayed. If \code{AllPlot==TRUE}, the x/y coordinates will be plotted as histograms on the x and y axis of the plot.}
\item{$histograms}{1D density plots of the spot localization on the length axis of the cells, grouped by \code{bycol} in \code{x==inp} groups.}
\item{$spotdata}{data frame containing the groups, relative x/y localizations of the spots and cell length/width}
\item{$meshdata}{if \code{MESH} is included, the function returns a dataframe including the rotated meshes and groups.}
}

\references{

%% ~put references to the literature/web site here ~
}
\author{
Renske van Raaphorst
}
\note{
This function is meant to give an overview of the data. To modify each individual plot, look at the shinyspots individual plot functions.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##use data from example dataset

##create list of plots, ordered by length, grouped in 6 groups.
GFPlist <- createPlotlist(GFP, 6, bycol="length", MESH, colorpalette="GreenYellow", AllPlot=T)

#have a look at the "plottotal" plot
plot(GFPlist$plottotal)

#save PDF of the 6 grouped plots (qplots)
ggplot2::ggsave(gridExtra::arrangeGrob(GFPlist$qplots, ncol=1, filename="GroupPlots.PDF"))

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
