\name{extr.MicrobeJ}
\alias{extr.MicrobeJ}
\alias{extr.MicrobeJMESH}
\alias{extr.MicrobeJSpots}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{upload MicrobeJ output into R for further analysis
%%  ~~function to do ... ~~
}
\description{
Upload the output of MicrobeJ into R. Upto now the code supports fluorescent localization data & cell outlines, both saved as seperate .CSV or .TXT files in the MicrobeJ GUI by first pressing \emph{<right click> -> expand all} on one of the Results columns and then pressing \emph{<right click> -> export} for both the cell dataset and the fluorescent localization dataset.
}
\usage{
extr.MicrobeJ(dataloc, spotloc)
}


\arguments{
  \item{dataloc}{path of the exported microbeJ comma-seperated file (.csv or .txt) containing the outlines of the bacterial cells.}
  \item{spotloc}{path of the exported microbeJ comma-seperated file (.csv or .txt) containing the sub-directory with fluorescent spot coordinates}

}
}
\details{
Input either \code{dataloc}, \code{spotloc} or both to receive a list containing one or two dataframes containing the microbeJ data.

}
\value{
Outlist containg:
\item{$MESH}{dataframe containing MicrobeJ outlines of the cell}
\item{$SPOTS}{dataframe containing MicrobeJ spot x/y coordinates and cell identifying number.}
}
\references{
\url{www.microbej.org}

Ducret, A., Quardokus, E.M. and Brun, Y.V., 2016. MicrobeJ, a tool for high throughput bacterial cell detection and quantitative analysis. \emph{Nature microbiology}, 1, p.16077.

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

##upload MicrobeJ datasets
dataloc <- file.choose()
spotloc <- file.choose()
##extract data
micdata <- extr.MicrobeJ(dataloc, spotloc)


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
