\name{plotTracks}
\alias{plotTracks}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot tracked spots and/or objects over time in 1 or 2 dimensions.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
plotTracks( meshdata,
            spotdata,
            objectdata,
            tracks = TRUE,
            ignore_singles = FALSE,
            movie = "none",
            timepalette_lines = "viridis",
            timepalette_fill = "magma",
            dimension = c("length", "width"),
            turn_cells = TRUE,
            mag,
            cell = "all",
            transparency = 0.2
            )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{meshdata}{
  A \code{\link{mesh}} dataframe containing cell outlines.
}
  \item{spotdata}{
  A \code{\link{Spot Data}} dataframe containing spot localizations.
}
  \item{objectdata}{
  A \code{\link{object}} dataframe containing fluorescent object outlines.
}
  \item{tracks}{
  When \code{tracks} == \code{TRUE}, the track of individual spots over time will be shown. When \code{tracks} == \code{FALSE}, only the spots will be shown.
}
  \item{ignore_singles}{
  When \code{ignore_singles} == \code{FALSE}, all detected spots will be shown. When \code{ignore_singles} ==\code{TRUE}, only those spots will be shown which are tracked over time.
}
  \item{movie}{
  When \code{movie}==\code{"none"} (default), the output will be a static plot. When \code{movie} == \code{"gganimate"},  \code{\href{https://github.com/thomasp85/gganimate}{gganimate}} (version => 1.0.3) will be used to make a GIF animation. Finally, when \code{movie} = \code{"plotly"}, \code{\href{https://plot.ly/}{plotly}} will be used to make an interactive HTML-based animation.
}
  \item{timepalette_lines}{
  A \code{\href{https://cran.r-project.org/web/packages/viridis/vignettes/intro-to-viridis.html}{viridis}} palette for both the spots & tracks as the outlines of the meshes, where color indicates the time in frames. Default = \code{"viridis"}.
}
  \item{timepalette_fill}{
A \code{\href{https://cran.r-project.org/web/packages/viridis/vignettes/intro-to-viridis.html}{viridis}} palette for the fill of the objects, where color indicates the time in frames. Default = \code{"magma"}.
}
  \item{dimension}{

}
  \item{turn_cells}{
%%     ~~Describe \code{turn_cells} here~~
}
  \item{mag}{
%%     ~~Describe \code{mag} here~~
}
  \item{cell}{
%%     ~~Describe \code{cell} here~~
}
  \item{transparency}{
%%     ~~Describe \code{transparency} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (meshdata, spotdata, objectdata, tracks = TRUE, ignore_singles = FALSE,
    movie = FALSE, timepalette_lines = "viridis", timepalette_fill = "magma",
    dimension = c("length", "width"), turn_cells = TRUE, mag,
    cell = "all", transparency = 0.2)
{
    outplot <- ggplot2::ggplot() + ggplot2::theme_minimal()
    if (cell[1] != "all" & is.numeric(cell) == FALSE) {
        stop("Cell number not recognized. Give cell number or keep default (cell='all') to return all cells")
    }
    if (missing(meshdata) != TRUE) {
        if ("length" \%in\% dimension & "width" \%in\% dimension) {
            if (cell[1] != "all") {
                meshdata <- meshdata[meshdata$cell \%in\% cell,
                  ]
            }
            if (turn_cells == TRUE) {
                if ("X_rot" \%in\% colnames(meshdata)) {
                  if ("Xrotum" \%in\% colnames(meshdata) != T) {
                    if (missing(mag)) {
                      stop("Please specify the pixel-micron conversion value 'mag'.")
                    }
                    meshdata$Xrotum <- meshdata$X_rot * unlist(get(magnificationList,
                      envir = magEnv)[mag])
                    meshdata$Yrotum <- meshdata$Y_rot * unlist(get(magnificationList,
                      envir = magEnv)[mag])
                  }
                }
                outplot <- outplot + ggplot2::geom_polygon(data = meshdata,
                  ggplot2::aes(x = Xrotum, y = Yrotum, color = frame,
                    group = frame), fill = "black", alpha = transparency/5)
            }
            if (turn_cells == FALSE) {
                outplot <- outplot + ggplot2::geom_polygon(data = meshdata,
                  ggplot2::aes(x = X, y = Y, color = frame, group = frame),
                  fill = "black", alpha = transparency/5)
            }
        }
    }
    if (missing(objectdata) != TRUE) {
        if (cell[1] != "all") {
            objectdata <- objectdata[objectdata$cell \%in\% cell,
                ]
        }
        if (turn_cells == TRUE) {
            if ("ob_out_x" \%in\% colnames(objectdata) != TRUE) {
                if ("ob_x" \%in\% colnames(objectdata) != TRUE) {
                  stop("Cannot find object coordinate data (ob_out_x/ob_out_y for turned cells, ob_x/ob_y for raw coordinates.")
                  if (missing(meshdata)) {
                    stop("Cannot find turned object data or mesh information to connect the objects to the cell localizations.\n                     Replace object dataframe for a dataframe (object_relative) containing this information or add mesh data to convert dataframe.")
                  }
                  if (missing(mag)) {
                    stop("Need pixel to micron conversion factor 'mag' to correctly connect mesh to object data.")
                  }
                  objectdata <- suppressWarnings(centrefun(objectframe))
                  objectdata <- suppressWarnings(midobject(meshdata,
                    objectframe, get(magnificationList, envir = magEnv)[mag]))
                }
            }
            outplot <- outplot + ggplot2::geom_polygon(data = objectdata,
                ggplot2::aes(x = ob_out_x, y = ob_out_y, group = (paste(obID,
                  frame)), fill = frame), color = "black", alpha = transparency)
        }
        if (turn_cells != TRUE) {
            if ("ob_x" \%in\% colnames(objectdata) != TRUE) {
                stop("Cannot find object coordinates 'ob_x'/'ob_y'.")
            }
            outplot <- outplot + ggplot2::geom_polygon(data = objectdata,
                ggplot2::aes(x = ob_x, y = ob_y, group = paste(obID,
                  frame), fill = frame), color = "black", alpha = transparency)
        }
        outplot <- outplot + scale_fill_viridis_c(option = timepalette_fill)
    }
    if (missing(spotdata) != TRUE) {
        if (cell[1] != "all") {
            spotdata <- spotdata[spotdata$cell \%in\% cell, ]
        }
        if (tracks == TRUE) {
            if (ignore_singles == TRUE) {
                spotdata <- spotdata[spotdata$trajectory != -1,
                  ]
            }
        }
        if (turn_cells == TRUE) {
            if ("Lmid" \%in\% colnames(spotdata) != T) {
                if ("l" \%in\% colnames(spotdata) != T) {
                  stop("Data doesn't include relative spot positions. Please use spotsInBox() to relate spot positions to cells.")
                }
                if (missing(mag)) {
                  stop("Please specify the pixel-micron conversion value 'mag'.")
                }
                spotdata$Lmid <- spotdata$l * unlist(get(magnificationList,
                  envir = magEnv)[mag])
                spotdata$Dum <- spotdata$d * unlist(get(magnificationList,
                  envir = magEnv)[mag])
            }
            outplot <- outplot + ggplot2::geom_point(data = spotdata,
                ggplot2::aes(x = Lmid, y = Dum, color = frame),
                size = 2)
            if (tracks == TRUE) {
                outplot <- outplot + ggplot2::geom_path(data = spotdata[spotdata$trajectory !=
                  -1, ], ggplot2::aes(x = Lmid, y = Dum, group = trajectory,
                  color = frame), size = 1)
            }
        }
        if (turn_cells == FALSE) {
            outplot <- outplot + ggplot2::geom_point(data = spotdata,
                ggplot2::aes(x = x, y = y, color = frame), size = 2)
            if (tracks == TRUE) {
                outplot <- outplot + ggplot2::geom_path(data = spotdata[spotdata$trajectory !=
                  -1, ], ggplot2::aes(x = x, y = y, group = trajectory,
                  color = frame), size = 1)
            }
        }
    }
    if (missing(meshdata) != T | missing(spotdata) != T) {
        outplot <- outplot + ggplot2::scale_color_viridis_c(option = timepalette_lines)
    }
    if (turn_cells == TRUE) {
        outplot <- outplot + ggplot2::coord_fixed() + facet_wrap(~cell)
    }
    if (turn_cells == FALSE) {
        outplot <- outplot + facet_wrap(~cell, scales = "free")
    }
    return(outplot)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
