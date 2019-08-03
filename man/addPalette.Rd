\name{addPalette}
\alias{addPalette}
\alias{showCurrentPalettes}
\alias{getPalette}
\alias{pal2Default}

\title{Manually add a heatmap palette and view the currently loaded palettes.

}
\description{\code{addPalette()}: function to manually add a colorpalette for heatmaps to the six heatmap palettes existing inside bactMAP. Use \code{showCurrentPalettes()} to print mockup-heatmap plots with the six default palettes and - if appliccable - the added palettes. Use \code{getPalette(palName)} to get the hex codes for the pallete you use. \code{pal2Default()} sets the color palettes back to default.
}
\usage{
addPalette(palList, palName)

showCurrentPalettes()

getPalette(palName)

pal2Default()

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{palList}{a list consisting of 3 hex color codes (character strings, for instance "#FFFFFF" or "#DD0088") or color names known by R (for instance "red" or "violet"). Optional; if argument is not used, the function will ask for the colors seperately while running.
}
 \item{palName}{a character string defining the name of the added colorpalette. Optional, if argument not used, the function will ask for the palette name.}
}

\value{\code{addPalette()} and \code{showCurrentPalettes()} will print heatmap plots using the color palettes loaded into bactMAP. \code{addPalette()} also adds the manually entered color palette to the list of color palettes loaded into the bactMAP environment during the R session.

}
\references{
The default color palettes are chosen based on the Color Universal Design (CUD, Masataka Okabe & Kei Ito): \url{http://jfly.iam.u-tokyo.ac.jp/color/}
}
\author{
Renske van Raaphorst
}

\seealso{
For more information and vast choice of good color palettes for plots, check \url{http://colorbrewer2.org}

For a comprehensive list of R color names, check \url{http://www.stat.columbia.edu/~tzheng/files/Rcolor.pdf}
}
\examples{
##Check the color palettes inside bactMAP
showCurrentPalettes()

##Pick three colors you like and put them in a list -
#in the order of lowest value - midpoint - highest value

listofcolors <- list("orangered", "darkseagreen", "burlywood1")

##Add to the color palettes inside bactMAP
addPalette(listofcolors, "MyBeautifulPalette")

##Get palette hex codes for "GreenYellow"
getPalette("GreenYellow")

##Put color palettes back to default
pal2Default()

}

\keyword{~color}
