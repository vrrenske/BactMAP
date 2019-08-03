\name{plotTreeBasic}
\alias{plotTreeBasic}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot cell genealogy information as networks using ggtree.
}
\description{
This function is using the ggtree package to show cell genealogy combined with cell information, e.g. fluorescence intensity. This part of bactMAP is under development - plotting cell trees is possible, but is not (yet) as straight-forward. Check the tutorials on \href{veeninglab.com/bactmap}{veeninglab.com/bactmap} to get more information on how to use this function.
}

\usage{
plotTreeBasic(phylo,
              extradata,
              yscalechange = FALSE,
              showClade = FALSE,
              layout = "rectangular",
              ydata,
              cellNumber,
              open.angle,
              linesize = 1,
              linecolor = "black",
              lines = TRUE,
              colors = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{phylo}{
  The \code{$generation_lists} output of \code{\link{extr_SuperSegger}} or \code{\link{extr_Oufti}}.
}
  \item{extradata}{
  Data connected to the \code{$generation_lists}: \code{$generation_dataframes}. When this is added, you can plot for instance the clades of your dataset. See the \href{veeninglab.com/bactmap}{tutorials} for more information.
}
  \item{yscalechange}{
  Default = \code{FALSE}. When \code{TRUE}, an y-axis is added to the tree. Only works when also indicating what needs to be plotted on the y-axis (indicate by putting the column name of choice from \code{extradata} as input in \code{ydata}).
}
  \item{showClade}{
  When \code{showClade==TRUE}, one clade will be highlighted. Decide which one by putting a cell number in \code{cellNumber}.
}
  \item{layout}{
  Inherited from \code{ggtree}; this decides the shape of the plotted tree. Options: \code{"rectangular", "slanted", "fan", "circular", "radial", "unrooted", "equal_angle", "daylight"}
}
  \item{ydata}{
  Indicate the column name of the variable you want to be plotted on the y-axis of the tree (when \code{yscalechange==TRUE}).
}
  \item{cellNumber}{
  Indicate the number of the cell of which you want to highlight the clade (when \code{showClade==TRUE}).
}
  \item{open.angle}{
  When \code{layout=="fan"}, indicate the size of the open angle (in degrees) here.
}
  \item{linesize}{
  Thickness of the lines in the tree.
}
  \item{linecolor}{
  Indicate the color of the lines in the tree.
  }
  \item{lines}{
  When highlighting a clade, set \code{lines} to \code{TRUE} when you want to do that by linetype.
}
  \item{colors}{
  When highlighting a clade, set \code{colors} to \code{TRUE} when you want to do that by color.
}
}

\value{
A \code{ggtree} object which can be edited using the \code{ggtree} commands.
}
\references{
  \href{https://yulab-smu.github.io/treedata-book/index.html}{The ggtree cookbook by Guangchuang Yu.}
}
\author{
Renske van Raaphorst
}

\examples{
\dontrun{
#this code requires other packages:
require(ggtree)
require(ggplot2)

# download example datasets from https://veeninglab.com/f/example_datasets.zip
load("ssbB_phylo.rda")
load("ssbB_meanfluo.rda")

#load example data from bactMAP package
myTreeData <- ssbB_phylo
myFluoData <- ssbB_meanfluo

#plot Tree with attributes
library(ggtree)

plotTreeBasic(myTreeData$phylos)\%<+\% myTreeData$data_attributes +
  geom_point(aes(x=x-branch.length, color=fluormean), size=6) +
  geom_tippoint(aes(color=fluormean_D), size=6) + #added tip points, different dot size
  geom_label(aes(x=x-branch.length, label=label), color=NA, fill="white", alpha=0.5, size=1) +
    #I also make a seethrough background label so you can see the text over the colored dots
  geom_text(aes(x=x-branch.length, label=label), size=3) +
  scale_color_viridis_c(option="magma") + #add different color scale
  theme(legend.position="right") + #add legend
  coord_flip()

}
}

