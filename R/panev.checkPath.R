#Script: panev.checkPath.R
#License: GPLv3 or later
#Modification date: 2019-05-13
#Written by: Valentino Palombo
#Contact: valentino.palombo@gmail.com
#Description: Check the pathway(s) of interest before run main PANEV scripts

panev.checkPath <- function(
  species=NULL,
  FL=NULL
)
{
  #check if your species code is correct
  panev.checkSpecies(species)
  #check the FL arguments
  if (is.null(FL)){
    cat("\n")
    cat("Pathway(s) is not specified. Please check! \n")
    stop("Exit",call. = F)
  }else{
    cat("\n")
    cat("Pathway(s) is specified... ")
  }
  #eliminate the species suffix from the pathway list
  orgs <- as.data.frame(KEGGREST::keggList("organism"))[2:3]
  KEGGpath <- data.table::setDT(as.data.frame(cbind(KEGGREST::keggList("pathway", species))), keep.rownames = T)
  string <- paste(" -", paste(orgs[which(orgs$organism == species),2]))    
  KEGGpath$V1 <- stringr::str_replace_all(KEGGpath$V1, stringr::fixed(string), "")
  #check if your pathway(s) of interest is present among the list of available pathways for your species
  for (i in 1:length(FL)){
    if (FL[i] %in% gsub(species, "map", KEGGpath$rn)){ 
      #cat("Your FL", paste(FL[i]), "pathway code is correct", sep=" ", "\n")
    }else{
      cat("\n")
      cat(paste(FL[i],"pathway code is incorrect or not available for your species. Please check!", sep=" ", "\n"))
      stop("Exit",call. = F)}
  }
  cat("and correct! ")
}