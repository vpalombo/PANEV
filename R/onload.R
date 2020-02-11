.onAttach <- function(libname, pkgname) {
  mymsg <- "Loading required package: PANEV\n\n\n"
  mymsg <- paste(mymsg,"Thanks for using PaNEV v1.0!\n")
  mymsg <- paste(mymsg,"For more information use: help(package = 'PANEV')\n")
  mymsg <- paste(mymsg,"                          browseVignettes(package = 'PANEV')\n\n")
  mymsg <- paste(mymsg,"Citation:\n")
  mymsg <- paste(mymsg,"  Authors: V. Palombo, M. Milanesi, G. Sferra, S. Capomaccio, S. Sgorlon, M. D'Andrea\n")
  mymsg <- paste(mymsg,"  Title: Pathways Network Visualization (PaNeV): an R package for a pathway-based network visualization\n")
  mymsg <- paste(mymsg,"  Journal: BMC Bioinformatics volume 21, Article number: 46 (2020)\n")
  mymsg <- paste(mymsg,"  DOI: https://doi.org/10.1186/s12859-020-3371-7\n")
  packageStartupMessage(mymsg)
  
  suppressWarnings(suppressPackageStartupMessages(library(bc3net)))
  suppressWarnings(suppressPackageStartupMessages(library(data.table)))
  suppressWarnings(suppressPackageStartupMessages(library(stringr)))
  suppressWarnings(suppressPackageStartupMessages(library(dplyr)))
  suppressWarnings(suppressPackageStartupMessages(library(XML)))
  suppressWarnings(suppressPackageStartupMessages(library(xml2)))
  suppressWarnings(suppressPackageStartupMessages(library(visNetwork)))
  suppressWarnings(suppressPackageStartupMessages(library(RColorBrewer)))
  suppressWarnings(suppressPackageStartupMessages(library(BiocManager)))
  suppressWarnings(suppressPackageStartupMessages(library(KEGGREST)))
  suppressWarnings(suppressPackageStartupMessages(library(biomaRt)))
  
  require(bc3net)
  require(data.table)
  require(stringr)
  require(dplyr)
  require(XML)
  require(xml2)
  require(visNetwork)
  require(RColorBrewer)
  require(BiocManager)
  require(KEGGREST)
  require(biomaRt)
  
}

