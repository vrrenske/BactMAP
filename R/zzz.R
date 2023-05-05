##load instructions


.onAttach <- function(libname,pkgname){
  packageStartupMessage("This is the development version of BactMAP updated on 5/5/2023. Use this version for the latest bug-fixes and new functions.
                        To download the stable version of bactmap, use 'remotes::install_github('veeninglab/bactmap')'")
}
