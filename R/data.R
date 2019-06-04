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

