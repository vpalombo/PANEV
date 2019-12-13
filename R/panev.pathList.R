#' @title PANEV pathways codes
#' @description The function retrieves the list of the available pathways in KEGG, according to a user defined search string. If path = \emph{NULL} the list of all available pathways is returned. To run all core PANEV functions a list of pathways of interest, coded with KEGG \emph{'path_ID'}, is required.
#' @usage panev.pathList(string = NULL)
#' @param string A string used to search within the description field of all KEGG available pathways (default = NULL).
#' @details This function is based on the \code{\link[KEGGREST]{keggList}} function of \pkg{KEGGREST} package (\url{http://bioconductor.org/packages/release/bioc/html/KEGGREST.html}.
#' @return A \emph{dataframe} containing two columns with KEGG \emph{path description} and \emph{path ID} respectively. 
#' @author Valentino Palombo (\email{valentino.palombo@gmail.com})
#' @references Tenenbaum D (2017). KEGGREST: Client-side REST access to KEGG. R package version 1.16.1. 
#' @examples ##### EXAMPLES CODE #####
#' #Create a list of all available pathways for PANEV
#' list <- panev.pathList(string = NULL)
#' 
#' #Look for a specific pathway(s) in PANEV, matching your search string
#' list <- panev.pathList(string = "lipid")


#Script: panev.pathList
#License: GPLv3 or later
#Modification date: 2019-05-13
#Written by: Valentino Palombo
#Contact: valentino.palombo@gmail.com
#Description: Get the list of pathways available for PANEV

panev.pathList <- function(
  string = NULL
)
{
  #download the KEGG pathway list
  path <- data.frame(as.matrix(KEGGREST::keggList("pathway", "")))
  path$path_ID <- row.names(path)
  colnames(path) <- c("path_description", "path_ID")
  rownames(path) <- NULL
  #check for a specific pathway(s) of interest
  if (is.null(string)){
    cat("The list of all available pathways was created! \n")
    cat("Remember to use the correct path Id(s) for relative PANEV functions. \n")
    return(path)
  }else{
    #look for the specified string
    name <- grep(string, path$path_description, ignore.case = T, value=T)
    path <- data.frame(path[path$path_description %in% name,], row.names = NULL)
    if (length(path$path_ID)>0){
      cat("The list of the pathway(s), matching your string, was created! \n")
      cat("Remember to use the correct path Id(s) for relative PANEV functions. \n")
      return(path)
    }else{
      cat("No match found for the requested string! \n")
    }
  }
}