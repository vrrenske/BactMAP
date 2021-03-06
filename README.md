BacTMAP
================

## BacTMAP - Bacteria Tool for Microscopy Analysis & Plotting

<img src="man//figures//logo_bactmap_pink.png" width="200" height="200" />

This package is meant to make it easier for microbiologists to combine
and analyse segmentation & fluorescence data derived from custom
software like Oufti, Morphometrics, MicrobeJ or Supersegger. There are
standard functions for importing several popular programs, but also
options for importing CSVs or Matlab files with a specific structure. If
your favorite segmentation or spot detection program is missing, don’t
hesitate to ask\!

I use ggplot and many ggplot-extensions to make visual summaries of the
segmentation & fluorescence data, show the cell genealogy, plot protein
trajectories in the cell and many other things. This can give you a
quick overview of your data to make a decision on further analysis and
custom visualization.

**This is the development page of BactMAP. Please go to https://github.com/veeninglab/bactmap for 
the stable build. There you can also find the wiki with instructions, documentation and tutorials.
If you have any comments, requests or find bugs, please report it on the [Issues page](https://github.com/veeninglab/bactmap/issues).**

## Download and install package

    #install devtools if not done yet
    install.packages("devtools")
    
    #install bactmap from my github repository
    devtools::install_github("vrrenske/bactMAP")
    
    #or from the stable build:
    devtools::install_github("veeninglab/bactMAP")
    
    #or install it with all dependencies
    devtools::install_github("vrrenske/bactMAP", dependencies=TRUE)
    
    #load package
    library(bactMAP)
