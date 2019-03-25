\name{Pixels2um}
\alias{Pixels2um}
\alias{Pixel2um}
\alias{pixels2um}
\alias{pixel2um}
\alias{getPixels2um}
\alias{addPixels2um}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Add and look up magnification factor "pixels2um"
%%  ~~function to do ... ~~
}
\description{Functions to check the pixel to micron conversion factors loaded into bactMAP and to add new conversion factors. The new conversion factors will be saved throughout the running R session.
}
\usage{
getPixels2um()

addPixels2um(pixelName, pixels2um)
}

\arguments{
  \item{pixelName}{A character string identifying your new conversion factor
  }
  \item{pixels2um}{A numerical value to convert the pixel values to microns}

  }

\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{getPixels2um() and addPixels2um("ConversionName", conversionfactor) both return the currently loaded conversion factors as a list of vectors, where the first item of each vector is the conversionfactor-name and the second item is the numeric conversion factor.
}
\references{
}
\author{
Renske van Raaphorst
}
\note{If you don't know your pixel conversion factor, open your original images in FIJI/ImageJ and check the Properties of the image.
}


\examples{

## To check the currently loaded pixel conversion factors:
getPixels2um()

## To add your pixel conversion factor to bactMAP:
addPixels2um("MyPixel2um", 0.543)

}
