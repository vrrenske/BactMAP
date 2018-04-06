---
title: "BacTMAP"
output: github_document
---



## BacTMAP - Bacteria Tool for Microscopy Analysis & Plotting

This R package is developed for easy processing, analysis & plotting of segmentation & fluorescence data from multiple popular bacterial (phase-contrast) segmentation tools, like Oufti, MicrobeJ & Morphometrics. The data can be combined with the raw data or spot detection data from other programs.

## Download and install package

```
#install devtools if not done yet
install.packages("devtools")

#install bactmap from my github repository
devtools::install_github("vrrenske/shinyspots")

#load package
library(bactMAP)
```

## Importing data

All importing functions start with **extr.**. See help(bactMAP::extr.all) in the package documentation for a full overview of all import functions.

In this case, the outlines of the cells and the ftsZ filaments where segmented using Oufti (www.oufti.org). BacTMAP has a function to directly upload oufti's matlab-file output into R using the **R.matlab**-package. However, in the case of object-detection, the matlab files are of a version incompatible with R.matlab. Therefore I saved my Oufti-output as CSV-files. 

BacTMAP's **extr.Oufti()** recognizes both .CSV and .MAT-files (under version 5.0). 

The spot localizations were fitted using the imageJ-plugin Peak Fitter from the package iSBatch (http://singlemolecule.github.io/iSBatch/). To upload the spot localizations, use **extr.ISBatch()**.



```r
#Add your own pixels to um conversion factor:
addPixels2um(pixels2um=0.0499538, pixelName="MyPixelConversion")
```

```
## [1] "Currently loaded magnification converters:"
## $`100x_LeicaVeening`
## [1] 0.0499538
## 
## $`100x_DVMolgen`
## [1] 0.06455
## 
## $No_PixelCorrection
## [1] 1
## 
## $MyPixelConversion
## [1] 0.0499538
## 
## $MyPixelConversion
## [1] 0.0499538
```

```r
#Load data into R
oufti_example <- system.file("extdata", "oufti_MK396_ftsz.csv", package = "bactMAP")
oufti_data <- extr.Oufti(oufti_example, "MyPixelConversion")
```

```
## [1] "Finished Oufti Extraction. Data list includes:"
##                 Length Class      Mode   
## cellList        14     data.frame list   
## mesh            19     data.frame list   
## objectframe      5     data.frame list   
## pixel2um         1     -none-     numeric
## object_relative 12     data.frame list
```

```r
#isbatch_example <- system.file("extdata", "isbatch_MK396_dnax.csv")
#spot_data <- extr.ISBatch(isbatch_example)
```

## Plotting data
This dataset contains the outlines of Streptococcus pneumoniae D39 cell and the spot localizations of the clamp-loader of the replication machinery DnaX linked to a superfolder-GFP (dataset from van Raaphorst, Kjos & Veening, PNAS, 2017). To get a quick overview of the localization, you can automatically plot the data in different ways.


```r
#Plot data for quick overview

oufti_plots <- createPlotlist(REP = oufti_data$object_relative, inp = 4, MESH = oufti_data$mesh, mag="MyPixelConversion")
```

```
## [1] "Calculating mean cell outlines.."
## [1] "Finished calculating mean cell outlines"
## [1] "Done plotting."
```

```r
print(summary(oufti_plots))
```

```
##            Length Class      Mode
## lengthplot  9     gg         list
## widthplot   9     gg         list
## qplots      4     -none-     list
## plottotal   4     gtable     list
## histograms  4     -none-     list
## spotdata   16     data.frame list
## meshdata   20     data.frame list
```

## Output
The output is a list of different plots which are useful for exploring the data. Let's have a look at them in more detail:

#Length/Width Kymographs
Here the cells are ordered by cell length, and the spot localization (density plot) and cell poles (white lines) are plotted. This gives an indication of the localization over the cell cycle when there is no time resolution (!note that time and cell length are often no direct proxys of each other, so interpret with caution!).

```r
oufti_plots$lengthplot
```

![plot of chunk length and width histograms](figure/length and width histograms-1.pdf)

```r
oufti_plots$widthplot
```

![plot of chunk length and width histograms](figure/length and width histograms-2.pdf)

#The average cell - totalplot
Here, the average cell shape (white dots) & spot localization is plotted. Side graphs show histograms of length/width spot localization.

```r
gridExtra::grid.arrange(oufti_plots$plottotal)
```

![plot of chunk totalplot](figure/totalplot-1.pdf)

#Localization per group - qplots & histograms
Here, the cells are divided in groups (the amount is indicated by variable *inp* in the function **createPlotlist()**). By default, they are divided by cell length. This can give more insight in the distinct localization or shape of subgroups of cells.
![plot of chunk qplots](figure/qplots-1.pdf)

