\name{colorsCUD}
\alias{colorsCUD}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Return CUD color palette with color titles
}
\description{
Colorblind friendly color palette (https://jfly.uni-koeln.de/color/) saved as a list of hex values.
}
\usage{
colorsCUD(colornames)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{colornames}{
  (Optional) A vector of color titles, for instance the different microscopy channels ("GFP", "RFP") to link each color of the CUD palette to.
}
}

\value{
A list of color names
}
\references{
https://jfly.uni-koeln.de/color/
}

\examples{
colorsCUD()
}

\keyword{ ~color}
