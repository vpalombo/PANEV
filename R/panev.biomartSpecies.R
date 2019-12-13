#' @title Gene list preparation species codes
#' @description The function retrieves the correct \emph{organism code} from \emph{ensembl.org} website. This code is mandatory to run \code{\link[PANEV]{panev.dataPreparation}} and \code{\link[PANEV]{panev.exprdataPreparation}}, useful to generate a properly formatted gene list for PANEV, switching up-to-date gene identifiers from \emph{ensembl} to \emph{entrez} annotation (or vice versa).  
#' @usage panev.biomartSpecies(string = NULL)
#' @param string A string used to search the organism code within the full list. If \emph{NULL}, a list with all available organisms is returned (default = NULL).
#' @details This function is based on \code{\link[biomaRt]{listDatasets}} function from \pkg{biomaRt} package (\url{http://bioconductor.org/packages/release/bioc/html/biomaRt.html}).
#' @return A \emph{dataframe} with the list of organisms available for gene list preparation. The \emph{dataframe} contains three columns with \emph{organism_code}, \emph{description} and \emph{version} respectively.
#' @author Valentino Palombo (\email{valentino.palombo@gmail.com})
#' @references Mapping identifiers for the integration of genomic datasets with the R/Bioconductor package biomaRt. Steffen Durinck, Paul T. Spellman, Ewan Birney and Wolfgang Huber, Nature Protocols 4, 1184-1191 (2009).
#' @references BioMart and Bioconductor: a powerful link between biological databases and microarray data analysis. Steffen Durinck, Yves Moreau, Arek Kasprzyk, Sean Davis, Bart De Moor, Alvis Brazma and Wolfgang Huber, Bioinformatics 21, 3439-3440 (2005).
#' @examples ##### EXAMPLES CODE #####
#' #Create a list of all available organisms
#' list <- panev.biomartSpecies(string = NULL)
#' 
#' #Look for a specific organism matching a search string 
#' list <- panev.biomartSpecies(string = "cow")


#Script: panev.biomartSpecies.R
#License: GPLv3 or later
#Modification date: 2019-05-13
#Written by: Valentino Palombo
#Contact: valentino.palombo@gmail.com
#Description: Find the organism code among the list of those available on biomaRt.

panev.biomartSpecies <- function(
  string = NULL
  )
{
  #download the biomaRt ensembl dataset and show all the list of organisms available
  mart <- biomaRt::useMart("ensembl")
  orgs <- biomaRt::listDatasets(mart)
  colnames(orgs) <- c("organism_code", "description", "version")
  
  #check the argument
  if (is.null(string)){
    cat("The list of all available species for gene id convertion and data preparation was created! \n")
    cat("Remember to use the correct organism code for relative PANEV functions. \n")
    return(orgs)
  }else{
    #download the biomaRt ensembl dataset and show only the organism(s) matching the string of interest
    org <- grep(string, orgs$description, ignore.case = T, value=T)
    if (length(org)==0){
      cat("No match found for the requested string! \n")
    }else{
      orgs <- data.frame(orgs[orgs$description %in% org,1:3], row.names = NULL)
      colnames(orgs) <- c("organism_code", "description", "version")
      cat("The list of available species matched your string was created! \n")
      cat("Remember to use the correct organism code for relative PANEV functions. \n")
      return(orgs)  
    }
  }  
}