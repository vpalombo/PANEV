.onAttach <- function(libname, pkgname) {
  mymsg <- "Loading required package: PANEV\n\n\n"
  mymsg <- paste(mymsg,"Thanks for using PaNEV v1.0!\n")
  mymsg <- paste(mymsg,"For more information use: help(package = 'PANEV')\n")
  mymsg <- paste(mymsg,"                          browseVignettes(package = 'PANEV')\n\n")
  mymsg <- paste(mymsg,"Citation:\n")
  mymsg <- paste(mymsg,"  Authors: V. Palombo, M. Milanesi, G. Sferra, S. Capomaccio, S. Sgorlon, M. D'Andrea\n")
  mymsg <- paste(mymsg,"  Title: Pathways Network Visualization (PaNeV): an R package for a pathway-based network visualization\n")
  mymsg <- paste(mymsg,"  Journal: ...\n")
  mymsg <- paste(mymsg,"  DOI: ...\n")
  packageStartupMessage(mymsg)
  
  suppressWarnings(suppressPackageStartupMessages(library(bc3net)))
  suppressWarnings(suppressPackageStartupMessages(library(biomaRt)))
  suppressWarnings(suppressPackageStartupMessages(library(KEGGREST)))
  suppressWarnings(suppressPackageStartupMessages(library(data.table)))
  suppressWarnings(suppressPackageStartupMessages(library(stringr)))
  suppressWarnings(suppressPackageStartupMessages(library(dplyr)))
  suppressWarnings(suppressPackageStartupMessages(library(XML)))
  suppressWarnings(suppressPackageStartupMessages(library(visNetwork)))
  suppressWarnings(suppressPackageStartupMessages(library(limma)))
  suppressWarnings(suppressPackageStartupMessages(library(RColorBrewer)))
  suppressWarnings(suppressPackageStartupMessages(library(bindrcpp)))
  
  require(bc3net)
  require(biomaRt)
  require(KEGGREST)
  require(data.table)
  require(stringr)
  require(dplyr)
  require(XML)
  require(visNetwork)
  require(limma)
  require(RColorBrewer)
  require(bindrcpp)
  
}

