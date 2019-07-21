#' DnaX and FtsZ timelapse of D39 S.pneumoniae. Cells are tracked using Oufti and converted to "mesh" dataframe using BactMAP::extr_Oufti
#'
#' @source van Raaphorst, Renske, Morten Kjos, and Jan-Willem Veening. "Chromosome segregation drives division site selection in Streptococcus pneumoniae." Proceedings of the National Academy of Sciences 114.29 (2017): E5959-E5968.
#' @format A data frame with columns:
#' \describe{
#'  \item{X}{x-coordinate of cell outline point (pixels)}
#'  \item{Y}{y-coordinate of cell outline point (pixels)}
#'  \item{cell}{cell identifier. unique per frame.}
#'  \item{frame}{frame identifier. in this case, an image was taken every 20 seconds, so 1 frame corresponds to 20s.}
#'  \item{num}{number indicating the order of cell outline points (useful for plotting polygons)}
#'  \item{max.width}{maximum cell width in pixels}
#'  \item{length}{the point of cell length from pole 0}
#'  \item{steplength}{the distance (over the length axis) of this cell outline point from the previous cell outline point}
#'  \item{max.length}{the maximum cell length of the given cell}
#'  \item{xy}{indicates whether it's the right or the left side of the cell (inherited from oufti)}
#'  \item{area}{cell area in pixels(2)}
#'  \item{angle}{the angle of the cells length axis to the horizontal line of the image.}
#'  \item{Xmid}{the mid-point (pixel x-coordinate) of the cell}
#'  \item{Ymid}{the mid-point (pixel y-coordinate) of the cell}
#'  \item{X_rot}{the x-coordinate of the cell outline point (in pixels) when the cell length axis is horizontal, with mid-cell at (0,0)}
#'  \item{Y_rot}{the y-coordinate of the cell outline point (in pixels) when the cell length axis is horizontal, with mid-cell at (0,0)}
#'  \item{max_um}{the maximum cell length in micron}
#'  \item{maxwum}{the maximum cell width in micron}
#'  \item{Xrotum}{the x-coordinate of the cell outline point (in micron) when the cell length axis is horizontal, with mid-cell at (0,0)}
#'  \item{Yrotum}{the y-coordinate of the cell outline point (in micron) when the cell lenth axis is horizontal, with mid-cell at (0,0)}
#' }
#' @examples
#' \dontrun{
#'  DnaX_mesh
#' }
"DnaX_mesh"


#' DnaX and FtsZ timelapse of D39 S.pneumoniae. DnaX spots are tracked using ISbatch and converted using extr_ISBatch
#'
#' @source van Raaphorst, Renske, Morten Kjos, and Jan-Willem Veening. "Chromosome segregation drives division site selection in Streptococcus pneumoniae." Proceedings of the National Academy of Sciences 114.29 (2017): E5959-E5968.
#' @format A data frame with columns:
#' \describe{
#'  \item{x}{x-coordinate of detected DnaX spot}
#'  \item{Y}{y-coordinate of detected DnaX spot}
#'  \item{cell}{cell identifier. unique per frame.}
#'  \item{frame}{frame identifier. in this case, an image was taken every 20 seconds, so 1 frame corresponds to 20s.}
#'  \item{max.width}{maximum cell width in pixels}
#'  \item{max.length}{the maximum cell length of the given cell}
#'  \item{l}{relative localization of DnaX on the length axis of the cell.}
#'  \item{d}{relative localization of DnaX on the width axis (diameter) of the cell}
#'  \item{trajectory}{identifier of the DnaX track}
#'  \item{trajectory_length}{the time (in frames) the trajectory could be followed}
#'  \item{displacement_sq}{the displacement of the molecule (in pixels, as calculated by ISBatch)}
#' }
#' @examples
#' \dontrun{
#'  DnaX_tracks
#' }
"DnaX_tracks"


#' DnaX and FtsZ timelapse of D39 S.pneumoniae. FtsZ spots are tracked using ISbatch and converted using extr_ISBatch
#'
#' @source van Raaphorst, Renske, Morten Kjos, and Jan-Willem Veening. "Chromosome segregation drives division site selection in Streptococcus pneumoniae." Proceedings of the National Academy of Sciences 114.29 (2017): E5959-E5968.
#' @format A data frame with columns:
#' \describe{
#'  \item{x}{x-coordinate of detected dnaX spot}
#'  \item{Y}{y-coordinate of detected dnaX spot}
#'  \item{cell}{cell identifier. unique per frame.}
#'  \item{frame}{frame identifier. in this case, an image was taken every 20 seconds, so 1 frame corresponds to 20s.}
#'  \item{max.width}{maximum cell width in pixels}
#'  \item{max.length}{the maximum cell length of the given cell}
#'  \item{l}{relative localization of FtsZ on the length axis of the cell.}
#'  \item{d}{relative localization of FtsZ on the width axis (diameter) of the cell}
#'  \item{trajectory}{identifier of the FtsZ track}
#'  \item{trajectory_length}{the time (in frames) the trajectory could be followed}
#'  \item{displacement_sq}{the displacement of the molecule (in pixels, as calculated by ISBatch)}
#' }
#' @examples
#' \dontrun{
#'  FtsZ_tracks
#' }
"FtsZ_tracks"



#' DnaX and FtsZ timelapse of D39 S.pneumoniae. This dataframe contains the fluorescence and mesh information of the FtsZ channel of cell number 4. The data is retreived by uploading the TIFF of FtsZ-RFP using extr_OriginalStacks and subsequently running extr_OriginalCells with this and the dataset bactMAP::DnaX_mesh. For disk space reasons, only cell 4 is saved here. Get the whole dataset on veeninglab.com/bactmap.
#'
#' @source van Raaphorst, Renske, Morten Kjos, and Jan-Willem Veening. "Chromosome segregation drives division site selection in Streptococcus pneumoniae." Proceedings of the National Academy of Sciences 114.29 (2017): E5959-E5968.
#' @format A data frame with columns:
#' \describe{
#'  \item{x}{x-coordinate of the pixel}
#'  \item{y}{y-coordinate of the pixel}
#'  \item{cell}{cell identifier. unique per frame.}
#'  \item{frame}{frame identifier. in this case, an image was taken every 20 seconds, so 1 frame corresponds to 20s.}
#'  \item{max.width}{maximum cell width in pixels}
#'  \item{max.length}{the maximum cell length of the given cell}
#'  \item{pointN}{number of specific pixel. unique per cell and frame}
#'  \item{values}{intensity of the pixel}
#'  \item{X_rot}{pixel x coordinate when the cell is turned so that the length axis is horizontal and mid-cell is at 0,0.}
#'  \item{Y_rot}{pixel y coordinate when the cell is turned so that the length axis is horizontal and mid-cell is at 0,0.}
#'  \item{xt}{the point marking one corner of the square pixel (x-coordinate)}
#'  \item{yt}{the point marking one corner of the square pixel (y-coordinate)}
#'  \item{area}{the cell area (pixel squared)}
#' }
#' @examples
#' \dontrun{
#'  TurnedCell4
#' }
#'
"TurnedCell4"



#' Dataset retrieved from movie previously published in Moreno-Gamez et al., 2017 \href{https://www.ncbi.nlm.nih.gov/pubmed/29021534}{pubMED}. Competence is induced with CSP and ssbB-GFP is expressed. The timeing between cells differs slightly. Cells are tracked using SuperSegger.
#'
#' @source Moreno-Gámez, Stefany, et al. "Quorum sensing integrates environmental cues, cell density and cell history to control bacterial competence." Nature communications 8.1 (2017): 854.
#' @format A data frame with columns:
#' \describe{
#'  \item{node}{the number of the node in the genealogy tree}
#'  \item{birth}{birthframe of this cell}
#'  \item{cell}{cell identifier. unique per frame.}
#'  \item{death}{death frame (or frame the cell divides) of this cell}
#'  \item{edgelength}{the time (in frames) of the cell division}
#'  \item{fluorsum}{the sum of the total fluorescence inside the cell. inherited from SuperSegger (see \code{\href{https://github.com/wiggins-lab/SuperSegger/wiki}{the SuperSegger Wiki}} for more information)}
#'  \item{fluormean}{the mean of the total fluorescence inside the cell. inherited from SuperSegger (see \code{\href{https://github.com/wiggins-lab/SuperSegger/wiki}{the SuperSegger Wiki}} for more information)}
#'  \item{fluorsum_D}{the sum of the fluorescence at cell division. inherited from SuperSegger (see \code{\href{https://github.com/wiggins-lab/SuperSegger/wiki}{the SuperSegger Wiki}} for more information)}
#'  \item{fluormean_D}{the mean of the total fluorescence at cell division. inherited from SuperSegger (see \code{\href{https://github.com/wiggins-lab/SuperSegger/wiki}{the SuperSegger Wiki}} for more information)}
#'  \item{parent}{the node of the parent of this cell}
#'  \item{child1}{the cell number of the daughter cell of this cell}
#'  \item{child2}{the cell number of the other daughter cell of this cell}
#'  \item{root}{the root of the genealogy tree. by default, this is 0.}
#'  \item{nodelabel}{the cell number corresponding to the node}
#'
#' }
#' @examples
#' \dontrun{
#'  ssbB_meanfluo
#' }
#'
"ssbB_meanfluo"



#' Dataset retrieved from movie previously published in Moreno-Gamez et al., 2017 \href{https://www.ncbi.nlm.nih.gov/pubmed/29021534}{pubMED}. Competence is induced with CSP and ssbB-GFP is expressed. The timeing between cells differs slightly. Cells are tracked using SuperSegger.
#'
#' @source Moreno-Gámez, Stefany, et al. "Quorum sensing integrates environmental cues, cell density and cell history to control bacterial competence." Nature communications 8.1 (2017): 854.
#' @format A phylo object containing the genealogy information of this timelapse movie.
#' \describe{
#' Phylogenetic tree with 103 tips and 95 internal nodes.
#' }
#' @examples
#' \dontrun{
#'  ssbB_phylos
#' }
#'
"ssbB_phylos"



#' Dataset retrieved from movie previously published in Moreno-Gamez et al., 2017 \href{https://www.ncbi.nlm.nih.gov/pubmed/29021534}{pubMED}. Competence is induced with CSP and ssbB-GFP is expressed. The timeing between cells differs slightly. Cells are tracked using SuperSegger.
#'
#' @source Moreno-Gámez, Stefany, et al. "Quorum sensing integrates environmental cues, cell density and cell history to control bacterial competence." Nature communications 8.1 (2017): 854.
#' @format A phylo object containing the genealogy information of this timelapse movie.
#' \describe{
#' iGRAPH network dataset containing the same information as \code{\link{ssbB_phylos}}
#' }
#' @examples
#' \dontrun{
#'  plot(ssbB_network)
#' }
#'
"ssbB_network"


#' VanFL dataset: objects WT
#'
#' @source van Raaphorst, Renske, Morten Kjos, and Jan-Willem Veening. "Chromosome segregation drives division site selection in Streptococcus pneumoniae." Proceedings of the National Academy of Sciences 114.29 (2017): E5959-E5968
#' @format an "object_relative" data frame (\code{\link{object}}) describing the shape and localization of fluorescent vancomycin in Wild Type S. pneumoniae cells.
#' \describe{
#' Cells were stained with fluorescent vancomycin and imaged as described in \href{https://www.pnas.org/content/early/2017/06/30/1620608114}{van Raaphorst, Kjos & Veening, 2017}. The cells and fluorescent objects where segmented using Oufti and imported in R using BactMAP's extr_Oufti function.
#' }
#' @examples
#' \dontrun{
#' View(VanFL_objWT)
#' }
#'
"VanFL_objWT"


#' VanFL dataset: objects mapZ mutant
#'
#' \describe{Cells were stained with fluorescent vancomycin and imaged as described in \href{https://www.pnas.org/content/early/2017/06/30/1620608114}{van Raaphorst, Kjos & Veening, 2017}. The cells and fluorescent objects where segmented using Oufti and imported in R using BactMAP's extr_Oufti function.}
#' @source van Raaphorst, Renske, Morten Kjos, and Jan-Willem Veening. "Chromosome segregation drives division site selection in Streptococcus pneumoniae." Proceedings of the National Academy of Sciences 114.29 (2017): E5959-E5968
#' @format an "object_relative" data frame (\code{\link{object}}) describing the shape and localization of fluorescent vancomycin in S. pneumoniae cells where the gene \code{mapZ} was replaced by a chloramphenicol resistance marker (\code{DmapZ::cmR}).
#' @examples
#' \dontrun{
#' View(VanFL_objDM)
#' }
#'
"VanFL_objDM"


#' VanFL dataset: cell outlines (mesh) WT
#'
#' \describe{Cells were stained with fluorescent vancomycin and imaged as described in \href{https://www.pnas.org/content/early/2017/06/30/1620608114}{van Raaphorst, Kjos & Veening, 2017}. The cells and fluorescent objects where segmented using Oufti and imported in R using BactMAP's extr_Oufti function.}
#' @source van Raaphorst, Renske, Morten Kjos, and Jan-Willem Veening. "Chromosome segregation drives division site selection in Streptococcus pneumoniae." Proceedings of the National Academy of Sciences 114.29 (2017): E5959-E5968
#' @format a "mesh" data frame (\code{\link{mesh}}) describing the shape and size of Wild Type S. pneumoniae cells.
#' @examples
#' \dontrun{
#'  View(VanFL_meshWT)
#'  }
#'
"VanFL_meshWT"



#' VanFL dataset: cell outlines (mesh) DM
#'
#' \describe{Cells were stained with fluorescent vancomycin and imaged as described in \href{https://www.pnas.org/content/early/2017/06/30/1620608114}{van Raaphorst, Kjos & Veening, 2017}. The cells and fluorescent objects where segmented using Oufti and imported in R using BactMAP's extr_Oufti function.}
#' @source van Raaphorst, Renske, Morten Kjos, and Jan-Willem Veening. "Chromosome segregation drives division site selection in Streptococcus pneumoniae." Proceedings of the National Academy of Sciences 114.29 (2017): E5959-E5968
#' @format a "mesh" data frame (\code{\link{mesh}}) describing the shape and localization of fluorescent vancomycin in S. pneumoniae cells where the gene \code{mapZ} was replaced by a chloramphenicol resistance marker (\code{DmapZ::cmR}).
#' @examples
#' \dontrun{
#'  View(VanFL_objWT)
#'  }
#'
"VanFL_meshDM"

