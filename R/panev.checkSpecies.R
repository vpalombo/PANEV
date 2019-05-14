#Script: panev.checkSpecies.R
#License: GPLv3 or later
#Modification date: 2019-05-13
#Written by: Valentino Palombo
#Contact: valentino.palombo@gmail.com
#Description: Check the species of interest before run main PANEV scripts

panev.checkSpecies <- function(
  species=NULL
)
{
  #check the argument
  if (is.null(species)){
    cat("\n")
    cat("Species code not specified. Please check! \n")
    stop("Exit",call. = F)
  }else{
    if (length(species) == 1){
      cat("\n")
      cat("Species code specified... ")
    }else{
      cat("\n")
      cat("ONLY one species code at a time is accept. Please check! \n")
      stop("Exit",call. = F)
    }
  }
  #download the KEGG organisms list and check your species code is present inside the list
  orgs <- as.data.frame(KEGGREST::keggList("organism"))[2]
  if (species %in% orgs$organism){
    cat("and correct! \n")    
  }else{
    cat("but incorrect. Please check! \n")
    cat("Remember you can use the related 'panev.speciesCode' function to discover all species available for PANEV.\n")
    stop("Exit",call. = F)
  }
}