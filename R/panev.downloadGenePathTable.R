#Script: panev.downloadGenePathTable.R
#License: GPLv3 or later
#Modification date: 2019-05-13
#Written by: Valentino Palombo
#Contact: valentino.palombo@gmail.com
#Description: Download KEGG pathways id and related genes table

panev.downloadGenePathTable <- function(
  species = NULL
)
{
  #download pathways id table (with related genes)
  KEGGgenes <<- data.frame(cbind(rn=names(KEGGREST::keggLink(paste(species),"pathway")), 
                                 V1=KEGGREST::keggLink(paste(species),"pathway")), stringsAsFactors = F, row.names = NULL)
}