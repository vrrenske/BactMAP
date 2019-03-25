\name{Spot Data}
\alias{Spot Data}
\alias{spots}
\alias{spotframe}
\alias{spots_relative}

\docType{data}
\title{
Spot Data: dataframes describing fluorescent spot coordinates and (optional) trajectories
}
\description{
When fluorescent spots are detected, the output "spotframe" (and when mesh data is included "spots_relative") is created. This dataset contains the x and y coordinates of the spots, the image frame in which the spots occurs, and if cell information is given, the relative localization of the cells. When spots where tracked over time, the spot trajectory, squared displacement and track length are given.
}

\format{
  spotframe contains:
  \describe{
    \item{\code{x}}{a numeric vector giving the x coordinate of the spot on the image frame}
    \item{\code{y}}{a numeric vector giving the y coordinate of teh spot on the image frame}
    \item{\code{frame}}{a numeric vector giving the number of the image frame}
  }
  when spots where tracked, the following columns are added:
  \describe{
    \item{\code{trajectory}}{a numeric vector giving an index for the spot trajectory. when spot is not part of a trajectory, this vector == -1}
    \item{\code{displacement_sq}}{a numeric vector giving the spots squared displacement (in pixels^2)}
    \item{\code{trajectory_length}}{a numeric vector giving the length (in frames) of the spot trajectory}
    }

  when cell information is giving, the spot data is listed in "spots_relative", which include the above and the following:
  \describe{
    \item{\code{l or Lmid}}{a numeric vector giving the relative localization over the length axis of the cell. l is given in pixels, Lmid in micron}
    \item{\code{d or Dum}}{a numeric vector giving the relative localization over the width axis of the cell. d is giving in pixels, Dum in micron}
    \item{\code{cell}}{the index of the cell the spot belongs to}
    \item{\code{max.width}}{a numeric vector giving the length of the cell}
    \item{\code{max.length}}{a numeric vector giving the width of the cell}
    \item{\code{pip}}{when pip == 1, the spot is found inside a cell}
    \item{\code{spot}}{a numeric vector indexing the spots per cell}
    \item{\code{totalspot}}{the total amount of spots per cell}
    \item{\code{cellframe}}{optional: the cell index and frame index pasted together with a "." as separator, to give a unique index to each cell}
  }
}

