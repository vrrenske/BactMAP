\name{extr_ObjectJ}
\alias{extr_ObjectJ}
\alias{chainTracer}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Upload and transform ObjectJ bounding box data.
}
\description{
ObjectJ detects a bounding box around each cell, which it uses to calculate spot localization, cell width and cell length along the medial axis. It is possible to save the bounding box points manually by defining them in the Columns section of the ObjectJ results. \href{https://github.com/vrrenske/bactMAPdata/blob/master/objectJ_settings.ojj}{Here} you can download an ObjectJ Settings file where the bounding box columns are already created.

}
\usage{
extr_ObjectJ( dataloc,
              mag = "No_PixelCorrection",
              boundingBoxX = c("X1", "X2", "X3", "X4", "X5", "X6", "X7", "X8", "X9", "X10", "X11"),
              boundingBoxY = c("Y1", "Y2", "Y3", "Y4", "Y5", "Y6", "Y7", "Y8", "Y9", "Y10", "Y11"),
              turn_meshes = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataloc}{
  The path to the saved .CSV or .TXT file from ObjectJ.
}
  \item{mag}{
  Magnification conversion factor (see \code{\link{pixels2um}})
}
  \item{boundingBoxX}{
  The column names of the 11 X-coordinates of the bounding box of the cells. Put these columns in your ObjectJ analysis manually or download the \href{https://github.com/vrrenske/bactMAPdata/blob/master/objectJ_settings.ojj}{ObjectJ Settings here}.
}
  \item{boundingBoxY}{
  The column names of the 11 Y-coordinates of the bounding box of the cells. Put these columns in your ObjectJ analysis manually or download the \href{https://github.com/vrrenske/bactMAPdata/blob/master/objectJ_settings.ojj}{ObjectJ Settings here}.
}
  \item{turn_meshes}{
  Default TRUE. When TRUE, the bounding boxes are all turned so the mid-point is at [0,0] while the length axis of the cell is parallel to the X-axis.
  }
}
\details{
 ObjectJ can already convert it's data to micron for you! Therefore, the default setting for "mag" is "no_PixelCorrection". Check what your settings are before converting the data using BactMAP.
}
\value{
 \item{$cellList}{
 The original CSV as dataframe.
 }
 \item{$mesh}{
 Dataframe containing the X and Y coordinates of the bounding box and cell dimensions. When \code{turn_meshes=TRUE}, it also contains the coordinates of the turned meshes, cell angle (to the x axis of the image) and mid-point coordinates.
 }
 \item{GFPframe}{
 Dataset contains self-defined (by the user) fluorescence information from the ObjectJ eg septal position or fluorescence intensity.
 }
 \item{$chainframe}{
 Dataset contains information on cell chains: chain ID, length, width and if specified, fluorescence information.
 }
}
\references{
https://sils.fnwi.uva.nl/bcb/objectj/

Vischer, Norbert OE, et al. "Cell age dependent concentration of Escherichia coli divisome proteins analyzed with ImageJ and ObjectJ." Frontiers in microbiology 6 (2015): 586.
}
\author{
Renske van Raaphorst
}
\note{
It is possible to define your own column names in ObjectJ - but (apart from the bounding box columns) bactMAP recognises the standard names only. See the [example objectJ settings file](https://github.com/vrrenske/bactMAPdata/blob/master/objectJ_settings.ojj) for the standard column names.
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
function (dataloc, mag = "No_PixelCorrection", boundingBoxX = c("X1",
    "X2", "X3", "X4", "X5", "X6", "X7", "X8", "X9", "X10", "X11"),
    boundingBoxY = c("Y1", "Y2", "Y3", "Y4", "Y5", "Y6", "Y7",
        "Y8", "Y9", "Y10", "Y11"))
{
    outlist <- list()
    oj <- read.table(dataloc, header = T, sep = "\t")
    celldats <- data.frame(Frame = unique(oj$Frame))
    celldats <- data.frame(Frame = celldats[order(celldats$Frame),
        ])
    celldats$f <- c(1:nrow(celldats))
    oj <- merge(celldats, oj)
    oj$frame <- oj$f
    oj$f <- NULL
    oj$Frame <- NULL
    oj_chain <- oj[!is.na(oj$ChainAxis), ]
    oj_chain$chainID <- oj_chain$n
    if (length(unique(oj_chain$GFPfluor) > 1)) {
        oj_chain <- oj_chain[, c("chainID", "ChainAxis", "ChainDia",
            "frame", "GFPfluor", "GFPMax", "GFPMean", "GFPMid")]
    }
    else {
        oj_chain <- oj_chain[, c("chainID", "ChainAxis", "ChainDia",
            "frame")]
    }
    outlist$cellList <- oj
    outlist$chainframe <- oj_chain
    oj <- oj[!is.na(oj$CellAxis), ]
    oj2 <- reshape2::melt(oj, measure.vars = boundingBoxX, value.name = "Xum",
        variable.name = "index_X")
    oj2 <- reshape2::melt(oj2, measure.vars = boundingBoxY, value.name = "Yum",
        variable.name = "index_Y")
    oj2$index_X <- lapply(oj2$index_X, function(x) strsplit(as.character(x),
        "X")[[1]][[2]])
    oj2$index_Y <- lapply(oj2$index_Y, function(x) strsplit(as.character(x),
        "Y")[[1]][[2]])
    oj2$index_X <- as.numeric(oj2$index_X)
    oj2$index_Y <- as.numeric(oj2$index_Y)
    oj2 <- oj2[oj2$index_X == oj2$index_Y, ]
    oj2$num <- oj2$index_X
    oj2$max.length <- oj2$CellAxis
    oj2$cell <- oj2$n
    oj2$chainID <- oj2$ParentID
    oj3 <- unique(oj2[, c("cell", "max.length", "chainID", "GFPfluor",
        "GFPMean", "GFPMax", "GFPMid", "allSepta", "manSepta",
        "manTrace", "Thr", "startX", "startY", "barX", "barY",
        "barW", "barH", "SeptaPos", "frame")])
    outlist$GFPframe <- oj3
    oj2 <- oj2[, c("cell", "frame", "max.length", "chainID",
        "Xum", "Yum", "num")]
    p2um <- as.numeric(get(magnificationList, envir = magEnv)[mag])
    oj2$X <- oj2$Xum/p2um
    oj2$Y <- oj2$Yum/p2um
    oj2 <- oj2[order(oj2$frame, oj2$cell, oj2$num), ]
    outlist$mesh <- oj2
    return(outlist)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
