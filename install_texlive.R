if (!requireNamespace("tinytex", quietly = TRUE)) {
  if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes", repos = "https://cran.rstudio.com/", quiet = TRUE)
    cat('devtools installed\n')
  }
  remotes::install_github(c('yihui/tinytex'), quiet = TRUE)
  tinytex::install_tinytex()
}
