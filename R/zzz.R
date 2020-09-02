##load instructions


.onAttach <- function(libname,pkgname){
  packageStartupMessage("This is the development version of BactMAP updated on 9/2/2020. Use this version for the latest bug-fixes and new functions.
                        To download the stable version of bactmap, use 'remotes::install_github('veeninglab/bactmap')'")
}
