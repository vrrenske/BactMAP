if (!requireNamespace("tinytex", quietly = TRUE)) {
  if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes", repos = "https://cran.rstudio.com/", quietly = TRUE)
    cat('devtools installed\n')
  }
  remotes::install_github(c('yihui/tinytex'), quietly = TRUE)
  tinytex::install_tinytex()
  tinytex::tlmgr_install()
  tinytex::tlmgr_update()
}
