\name{micron}
\alias{micron}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Get back the unicode character for the micro-sign + m (micron).
}
\description{
Use micron() in a character string, for instance when changin the axis labels of a plot, to get the unicode for the micro-sign + m.
}
\usage{
micron()
}

\value{
character string: 'u00b5m'
}

\author{
Renske van Raaphorst
}

\examples{

paste("This is the sign for micron:", bactMAP::micron(), sep=" ")

}
