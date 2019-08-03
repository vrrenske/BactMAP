\name{perc_Division}
\alias{perc_Division}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Filter timelapse data by growth speed and whether cells underwent a whole division cycle.

}
\description{
  The function takes timelapse dataframes including cell size, frame and cell ID and return a dataframe with information on growth speed, whether the cell underwent a full division (and how many divisions), and a binning of each cell division in percentages. Next to this, it gives a summary of the parameters in the input dataframe per binned division point. It also plots the individual cell growth over binned division, as well as the average cell growth.
}
\usage{
perc_Division(timelapse, av = TRUE, plotgrowth = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{timelapse}{
  A dataframe from timelapse data containing at least the variables \code{frame}, \code{cell} and \code{max.length}.
}
  \item{av}{
  When \code{av}==\code{TRUE}, the function returns a dataframe with the variables of the input dataframe averaged per binned division moment ("percentage_binned").
}
  \item{plotgrowth}{
  When \code{plotgrowth}==\code{TRUE}, the function returns a plot of the cell growth by binned division point, growth speed and nth division.
}
}
\details{
  This function estimates whether a division is a "full" division by a set of cutoffs. In some cases, the segmentation data gives this information already. When this is the case, it can be more straightforward to use this data instead of this function.
}
\value{

1. A dataframe \code{timelapse} containing the following variables on top of the variables already in the input dataframe:
\item{division}{the nth division in the timelapse this cell is making}
\item{min_frame}{the frame number where this division starts}
\item{max_frame}{the frame where the cell is dividing}
\item{av_length}{the mean cell length over the whole division}
\item{length_var}{the variance of the cell length during the whole division}
\item{division_time}{the time this division takes (in frames)}
\item{fulldivision}{indicates whether the division is a full division (1 = full)}
\item{coeff}{the coefficient of the linear growth function of the cell for this division - indicator of growth speed}
\item{growth}{the divisions are grouped in 4 groups: no growth ("none"), fast, medium and slow growth.}
\item{percentage}{if the moment of cell birth is 0 and cell division is 100, the percentage gives how far the cell is from dividing at this moment}
\item{percentage_binned}{percentage grouped in 10 groups, for averaging}

2. A dataframe \code{mean_by_percentage}, where the variables are averaged over \code{percentage_binned}.

3. A ggplot \code{plot_growth} showing the cell size over the percentage of division.

4. A ggplot \code{plot_avgrowth} showing the average cell size over the percentage of division.
}

\author{
Renske van Raaphorst
}

\examples{
\dontrun{
#this example won't have much meaning because the timelapse was quite short.
#but analysis is quick because the dataset is not very big:
#for a better example, download our DnaX dataset from veeninglab.com/bactmap

#Download FtsZ_tracks.rda from https://veeninglab.com/f/example_datasets.zip

load("FtsZ_tracks.rda")

percD <- perc_Division(FtsZ_tracks)

percD$plot_avgrowth
}
}
