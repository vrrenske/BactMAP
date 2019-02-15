\name{combineDataframes}
\alias{combineDataframes}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to combine data by condition (e.g. mutant vs wild type) and/or channel (e.g. GFP/RFP), meant for BactMAP imported data.
}
\description{
The function adds a column "condition" and/or "channel" to each dataframe in the given *listofdataframes*, after which it concatenates all dataframes. If not all columns in the dataframes match, it discards the columns which don't match. The function returns the concatenated dataframe "finalframe" and the original list of dataframes with the "condition" and/or "channel" column added to each dataframe.
}
\usage{
combineDataframes(listofdataframes, listofconditions, listofchannels)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{listofdataframes}{
list of dataframes of which the common columns need to be combined
}
  \item{listofconditions}{
list of conditions corresponding to each dataframe. Needs to have the same order and length as listofdataframes.
}
  \item{listofchannels}{
list of channels corresponding to each dataframe. Needs to have the same order and length as listofdataframes.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
\item{finalframe}{concatenated dataframe with columns "condition" and/or "channel" added}
\item{originaldata}{the input list of dataframes with columns "condition" and/or "channel" added}
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

GFPdata_WT <- extr.Oufti(file.choose())
RFPdata_WT <- extr.Oufti(file.choose())

RFPdata_mut  <- extr.Oufti(file.choose())
GFPdata_mut  <- extr.Oufti(file.choose())

combined <- combineDataframes(list(GFPdata_WT, RFPdata_WT, RFPdata_mut, GFPdata_mut),
                              listofconditions = list("WT", "WT", "mut", "mut"),
                              listofchannels=list("GFP", "RFP", "RFP", "GFP"))

}
