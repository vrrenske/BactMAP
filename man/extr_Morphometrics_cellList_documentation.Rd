\name{extr.Morphometrics.cellList}
\alias{extr.Morphometrics.cellList}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to extract the complete cellList out of a morphometrics output .MAT file.
}
\description{
Use this function when you are interested in other data than the output of \emph{extr.Morphometrics()} filters out of the morphometrics output.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
extr.Morphometrics.cellList(morphpath)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{morphpath}{The path of a .MAT file containing the output of a morphometrics analysis.
%%     ~~Describe \code{x} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Dataframe containing all values of the morphometrics output file. Values depend on settings of morphometrics.

%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
\url{https://simtk.org/projects/morphometrics}

Ursell, T., Lee, T.K., Shiomi, D., Shi, H., Tropini, C., Monds, R.D., Colavin, A., Billings, G., Bhaya-Grossman, I., Broxton, M. and Huang, B.E., 2017. Rapid, precise quantification of bacterial cellular dimensions across a genomic-scale knockout library. \emph{BMC biology}, 15(1), p.17.
%% ~put references to the literature/web site here ~
}
\author{
Renske van Raaphorst
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{extr.Morphometrics}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##Choose .MAT morphometrics output:
morphpath <- file.choose()
##Get complete cell segmentation information:
cellList <- extr.Morphometrics.cellList(morphpath)
##Save as excel file using "xlsx" package
xlsx::write.xslx(cellList, "morphometrics_output.xlsx", sheetName="Sheet1")
}

