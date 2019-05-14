#Script: panev.downloadPathDescription.R
#License: GPLv3 or later
#Modification date: 2019-05-13
#Written by: Valentino Palombo
#Contact: valentino.palombo@gmail.com
#Description: Download the list of available KEGG pathways with the related description

panev.downloadPathDescription <- function(
  species = NULL
)
{
  #download KEGG pathways description
  KEGGpath <<- data.table::setDT(as.data.frame(cbind(KEGGREST::keggList("pathway", species))), keep.rownames = T)
  orgs <- as.data.frame(KEGGREST::keggList("organism"))[2:3]
  #eliminate the species suffix from the list
  string <- paste(" -", paste(orgs[which(orgs$organism == species),2]))    
  KEGGpath$V1 <<- stringr::str_replace_all(KEGGpath$V1, stringr::fixed(string), "")
}