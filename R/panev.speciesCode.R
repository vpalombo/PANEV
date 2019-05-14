#' @title PANEV species codes
#' @description Get the list of species available for panev. The correct \emph{'panev code'} is mandatory for run correctly the main PANEV functions.
#' @usage panev.speciesCode(string = NULL)
#' @param string A string used to search the organism code within the full list. (default = NULL).
#' @details This function is based on the \code{\link[KEGGREST]{keggList}} function of \pkg{KEGGREST} package (\url{http://bioconductor.org/packages/release/bioc/html/KEGGREST.html}.
#' @return A \emph{dataframe} containing two columns named \emph{'species'} and \emph{'panev code'}, respectively.
#' @author Valentino Palombo (\email{valentino.palombo@gmail.com})
#' @references Tenenbaum D (2017). KEGGREST: Client-side REST access to KEGG. R package version 1.16.1. 
#' @examples ##### EXAMPLES CODE #####
#' #Create a list of all available species for PANEV
#' list <- panev.speciesCode()
#' 
#' #Look for a specific species in PANEV, matching a search string
#' list <- panev.speciesCode(string="bos")


#Script: panev.speciesCode
#License: GPLv3 or later
#Modification date: 2019-05-13
#Written by: Valentino Palombo
#Contact: valentino.palombo@gmail.com
#Description: Show the list of species available for PANEV

panev.speciesCode <- function(
  string = NULL
  )
{
  #download the list of KEGG available organisms
  orgs <- as.data.frame(KEGGREST::keggList("organism"))[2:3]
  #orgs <- as.data.frame(cbind(paste(orgs$organism), paste(orgs$species)))
  colnames(orgs) <- c("panev_code", "species" )
  #check for a specific species of interest
  if (is.null(string)){
    cat("The list of all species available for PANEV analysis was created! \n")
    cat("Remember to use the correct organism code for relative PANEV functions. \n")
    return(orgs)    
  }else{
    #look for the specified string
    species <- grep(string, orgs$species, ignore.case = T, value=T)
    if (length(species) >= 1){
      panev_code <- paste(orgs[orgs$species %in% species,1])
      result <- as.data.frame(cbind(species,panev_code))
      cat("The list of available species, matched your string, was created! \n")
      cat("Remember to use the correct organism code for relative PANEV functions. \n")
      return(result)  
    }else{
      cat("No match found for the requested string! \n")
    }
  }
}