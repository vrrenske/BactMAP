##load instructions


.onAttach <- function(libname,pkgname){
  packageStartupMessage("This is the development version of BactMAP. To download the stable version of bactmap, use 'remotes::install_github('veeninglab/bactmap')'")
}
