#' @title PANEV gene expression dataset preparation 
#' @description In case of an expression dataset, to run \code{\link[PANEV]{panev.exprscript}} function, PANEV requires a \emph{dataframe} containing five columns labelled as follows: \emph{ensembl gene id}, \emph{entrezgene}, \emph{external gene name}, \emph{FC}, \emph{pvalue}. 
#' This function helps to create a properly formatted dataframe, converting the gene list from \emph{ensembl} or \emph{entrez} gene annotation (or vice versa) and adding the \emph{gene symbol}.
#' @usage panev.exprdataPreparation(in.file, gene_id = NULL, biomart.species = NULL)
#' @param in.file Name of input file (with extension). The input file is a table with three columns labelled \emph{ensembl_gene_id} OR \emph{'entrezgene'}, \emph{'FC'} and \emph{'pvalue'}, respectivelly. The input file \bold{must} rely in the working directory.
#' @param gene_id Type of annotation of provided gene list. You can choose among \emph{ensembl} or \emph{entrez} gene ID annotation (default = NULL).
#' @param biomart.species The biomaRt species code of interest. You can check the list of all available species on biomaRt with the function \code{\link[PANEV]{panev.biomartSpecies}} (default = NULL).
#' @details This function is based on the main \pkg{biomaRt} (\url{http://bioconductor.org/packages/release/bioc/html/biomaRt.html} query function \code{\link[biomaRt]{getBM}} which, given a set of filters and corresponding values, it retrieves the user specified attributes from the BioMart database one is connected to.
#' @return A \emph{<in.file>_converted.txt} file with five columns labeled \emph{'ensembl_gene_id'}, \emph{'entrezgene'} and \emph{'external_gene_name'}, \emph{'FC'}, \emph{'pvalue'}, stored in the working directory. 
#' @author Valentino Palombo (\email{valentino.palombo@gmail.com})
#' @references Mapping identifiers for the integration of genomic datasets with the R/Bioconductor package biomaRt. Steffen Durinck, Paul T. Spellman, Ewan Birney and Wolfgang Huber, Nature Protocols 4, 1184-1191 (2009).
#' @references BioMart and Bioconductor: a powerful link between biological databases and microarray data analysis. Steffen Durinck, Yves Moreau, Arek Kasprzyk, Sean Davis, Bart De Moor, Alvis Brazma and Wolfgang Huber, Bioinformatics 21, 3439-3440 (2005).
#' @examples ##### EXAMPLES CODE #####
#' #Copy the example files in the current working directory
#' panev.example()
#' 
#' #Look for the organism code matching the search string 
#' list <- panev.biomartSpecies(string = "pig")
#' biomart.species <- as.character(list[4,1]) # sscrofa_gene_ensembl
#' 
#' #Convert the gene list with ensembl gene id
#' # Copy the example data file 'gene_list.txt' in the current working directory
#' genelist_converted <- panev.exprdataPreparation(in.file = "ensembl_expr_genelist.txt", 
#'                                               gene_id="ensembl", 
#'                                               biomart.species = biomart.species)
#' 
#' #Convert the gene list with entrez gene id
#' # Copy the example data file 'gene_list.txt' in the current working directory
#' genelist_converted <- panev.exprdataPreparation(in.file = "entrez_expr_genelist.txt", 
#'                                               gene_id="entrez", 
#'                                               biomart.species = biomart.species)
#' 

#Script: panev.exprdataPreparation.R
#License: GPLv3 or later
#Modification date: 2019-05-13
#Written by: Valentino Palombo
#Contact: valentino.palombo@gmail.com
#Description: Convert gene ID from ensembl to entrez (or vice versa) and to add the gene symbol to expression datasets.

panev.exprdataPreparation <- function(
  in.file,
  gene_id=NULL,
  biomart.species=NULL
)
{
  #check the arguments genelist, biomartspecies, gene_id
  if (file.exists(in.file)){
    genelist <- read.table(in.file, header=T, stringsAsFactors = F)
    cat("Input file imported! \n")
  }else{
    cat("Input file not found. Please check! \n")
    stop("Exit",call. = F)
  }
  
  if (is.null(biomart.species)){
    cat("\n")
    cat("BiomaRt species not specified. Please check! \n")
    stop("Exit",call. = F)
  }else{
    orgs <- biomaRt::listDatasets(biomaRt::useMart("ensembl"))
    if (biomart.species %in% orgs$dataset){
      cat("\n")
      cat("BiomaRt species correct! \n")
    }else{
      cat("\n")
      cat("BiomaRt species incorrect. Please check! \n")
      stop("Exit",call. = F)
    }
  }
  
  if (is.null(gene_id)){
    cat("\n")
    cat("Gene id not specified. Please check! \n")
    stop("Exit",call. = F)
  }else{    
    if (gene_id %in% c("ensembl","entrez")){
      if ((gene_id == "ensembl" & ("ensembl_gene_id" %in% colnames(genelist))) | 
          (gene_id == "entrez" & ("entrezgene" %in% colnames(genelist)))){
            cat("\n")
            cat("Gene id correct! \n")
          }else{
            cat("\n")
            cat("Gene id does not macth with column label. Please check! \n")
            stop("Exit",call. = F)
          }
    }else{
      cat("\n")
      cat("Gene id incorrect. You must choose among <ensembl> or <entrez>. \n")
      stop("Exit",call. = F)
    }
  }
  #generate output file name
  in.file <- gsub("\\..*", "", in.file)
  #convertion gene ID
  orgs <- biomaRt::listDatasets(biomaRt::useMart("ensembl"))
  if (gene_id=="ensembl"){
    cat("\n")
    cat("Converting from ensembl ID to entrez ID ... \n")
    converted <- biomaRt::getBM(attributes = c("ensembl_gene_id", "entrezgene", "external_gene_name"),
                                filters = c("ensembl_gene_id"),
                                values = genelist,
                                mart=(biomaRt::useMart("ensembl", biomart.species)), verbose = FALSE)
    genelist <- merge(genelist, converted, by.x = "ensembl_gene_id", by.y = "ensembl_gene_id", all = F)
    genelist <- unique(genelist)
    cat("DONE\n")
    #cat("\n")
    #cat("Checking for corresponding genes in KEGG database ... \n")
    cat("n.", paste(sum(length(which(!is.na(genelist$entrezgene))))), "of your", paste(sum(length(which(!is.na(genelist$ensembl_gene_id))))), "genes have corresponding gene in KEGG database. \n", sep= " ")
    write.table(x = genelist, file = paste(in.file,"_converted.txt",sep=""), row.names = F)
    cat("\n")
    cat("Gene list exported!")
    return(genelist)
  }else{
    cat("\n")
    cat("Converting from entrez ID to ensembl ID ... \n")
    converted <- biomaRt::getBM(attributes = c("ensembl_gene_id", "entrezgene", "external_gene_name"),
                                filters = c("entrezgene"),
                                values = genelist,
                                mart=(biomaRt::useMart("ensembl", biomart.species)), verbose = FALSE)
    genelist <- merge(genelist, converted, by.x = "entrezgene", by.y = "entrezgene", all = F)
    genelist <- unique(genelist)
    cat("DONE\n")
    #cat("\n")
    #cat("Checking for corresponding genes in KEGG database ... \n")
    cat("n.", paste(sum(length(which(!is.na(genelist$entrezgene))))), "of your", paste(sum(length(which(!is.na(genelist$ensembl_gene_id))))), "genes have corresponding gene in KEGG database. \n", sep= " ")
    write.table(x = genelist, file = paste(in.file,"_converted.txt",sep=""), row.names = F)
    cat("\n")
    cat("Gene list exported!")
    return(genelist)
  }
}
