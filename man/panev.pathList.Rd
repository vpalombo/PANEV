% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/panev.pathList.R
\name{panev.pathList}
\alias{panev.pathList}
\title{PANEV pathways codes}
\usage{
panev.pathList(string = NULL)
}
\arguments{
\item{string}{A string used to search within the description field of all KEGG available pathways (default = NULL).}
}
\value{
A \emph{dataframe} containing two columns with KEGG \emph{path description} and \emph{path ID} respectively.
}
\description{
The function retrieves the list of the available pathways in KEGG, according to a user defined search string. If path = \emph{NULL} the list of all available pathways is returned. To run all core PANEV functions a list of pathways of interest, coded with KEGG \emph{'path_ID'}, is required.
}
\details{
This function is based on the \code{\link[KEGGREST]{keggList}} function of \pkg{KEGGREST} package (\url{http://bioconductor.org/packages/release/bioc/html/KEGGREST.html}.
}
\examples{
##### EXAMPLES CODE #####
#Create a list of all available pathways for PANEV
list <- panev.pathList(string = NULL)

#Look for a specific pathway(s) in PANEV, matching your search string
list <- panev.pathList(string = "lipid")
}
\references{
Tenenbaum D (2017). KEGGREST: Client-side REST access to KEGG. R package version 1.16.1.
}
\author{
Valentino Palombo (\email{valentino.palombo@gmail.com})
}
