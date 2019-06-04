\name{plotCellsTime}
\alias{plotCellsTime}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Function to plot cell fluorescence from a TIFF over time.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
plotCellsTime(celdat, updown = T, movie = F, viridisoption = "magma", cellN, minf, maxf)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{celdat}{
  TIFF and mesh outline data from a timelapse movie as it was obtained by running \code{extr_OriginalCells}.
}
  \item{updown}{
  When \code{updown}==\code{TRUE}, a horizontal cell tower is made, otherwise, the cells are displayed left-right. \code{movie}=={TRUE} will overwrite this and return an animation instead.
}
  \item{movie}{
When \code{movie} == \code{TRUE}, an animation made with \code{\href{https://github.com/thomasp85/gganimate}{gganimate}} is returned.
}
  \item{viridisoption}{
  The color palette used. The default palettes are part of the \code{\href{https://cran.r-project.org/web/packages/viridis/vignettes/intro-to-viridis.html}{viridis}} color palettes (default = "magma").
}
  \item{cellN}{
  The cell number or numbers to be displayed. When not indicated, a list of plots of all cells will be made.
}
  \item{minf}{
  The minimum frame to be displayed.
}
  \item{maxf}{
  The maximum frame to be displayed
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A plot, list of plots or animation of a single cell's fluorescence over time.
}
\references{
%% ~put references to the literature/web site here ~
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

##In BactMAP, the raw data of one cell is saved as "TurnedCell4":

#to plot a left-right tower:
plotCellsTime(bactMAP::TurnedCell4, updown=F, viridisoption="inferno")

#or to create an animation:
plotCellsTime(bactMAP::TurnedCell4, movie=T)

}
