#' @title PANEV gene list preparation 
#' @description The function helps to create a properly formatted gene list required by \code{\link[PANEV]{panev.network}}: a \emph{dataframe} with three columns labelled as '\emph{ensembl_gene_id}', '\emph{entrezgene}' and '\emph{external_gene_name}', respectively.
#' @usage panev.dataPreparation(in.file, gene_id = NULL, biomart.species = NULL)
#' @param in.file Name of input file (with extension). The input file is a one column list of genes of interest, in \emph{ensembl} or \emph{entrez} gene annotation. The file \bold{must} rely in the working directory.
#' @param gene_id Type of gene identifiers provided. User can choose among \emph{ensembl} or \emph{entrez gene ID} annotation (default = NULL).
#' @param biomart.species The biomaRt organism code. The handy function \code{\link[PANEV]{panev.biomartSpecies}} is provided to retrieve the correct code (default = NULL).
#' @details This function is based on the main \pkg{biomaRt} (\url{http://bioconductor.org/packages/release/bioc/html/biomaRt.html}) query function \code{\link[biomaRt]{getBM}} which, given a set of filters and corresponding values, it retrieves the user specified attributes from the BioMart database one is connected to.
#' @return A \emph{<in.file>_converted.txt} file containing three columns with \emph{'ensembl_gene_id'}, \emph{'entrezgene'} and \emph{'external_gene_name'} respectively, stored in the working directory. 
#' @author Valentino Palombo (\email{valentino.palombo@gmail.com})
#' @references Mapping identifiers for the integration of genomic datasets with the R/Bioconductor package biomaRt. Steffen Durinck, Paul T. Spellman, Ewan Birney and Wolfgang Huber, Nature Protocols 4, 1184-1191 (2009).
#' @references BioMart and Bioconductor: a powerful link between biological databases and microarray data analysis. Steffen Durinck, Yves Moreau, Arek Kasprzyk, Sean Davis, Bart De Moor, Alvis Brazma and Wolfgang Huber, Bioinformatics 21, 3439-3440 (2005).
#' @examples ##### EXAMPLES CODE #####
#' #Copy the example files in the current working directory
#' panev.example()
#' 
#' #Look for the organism code matching the search string 
#' list <- panev.biomartSpecies(string = "cow")
#' biomart.species <- as.character(list[1,1]) # btaurus_gene_ensembl
#' 
#' #Prepare PANEV input file using a gene list containing ensembl gene id.
#' genelist_converted <- panev.dataPreparation(in.file = "ensembl_genelist.txt", 
#'                                           gene_id = "ensembl", 
#'                                           biomart.species = biomart.species)
#' 
#' #Prepare PANEV input file using a gene list containing entraz gene id.
#' genelist_converted <- panev.dataPreparation(in.file = "entrez_genelist.txt", 
#'                                           gene_id = "entrez", 
#'                                           biomart.species = biomart.species)
#' 

#Script: panev.dataPreparation.R
#License: GPLv3 or later
#Modification date: 2019-05-13
#Written by: Valentino Palombo
#Contact: valentino.palombo@gmail.com
#Description: Prepare the dataset for PANEV core function. This function converts gene ID from ensembl to entrez (or vice versa) and to add the gene symbol.

panev.dataPreparation <- function(
  in.file,
  gene_id=NULL,
  biomart.species=NULL
  )
{
  #check the arguments genelist, biomartspecies, gene_id
  if (file.exists(in.file)){
    genelist <- read.table(in.file, header=F, stringsAsFactors = F)
    cat("Input file imported! \n")
  }else{
    cat("Input file not found. Please check! \n")
    stop("Exit",call. = F)
  }
  
  if (is.null(biomart.species)){
    cat("\n")
    cat("Biomart species not specified. Please check! \n")
    stop("Exit",call. = F)
  }else{
    orgs <- biomaRt::listDatasets(useMart("ensembl"))
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
      cat("\n")
      cat("Gene id correct! \n")
    }else{
      cat("\n")
      cat("Gene id incorrect. You must choose among <ensembl> or <entrez>. \n")
      stop("Exit",call. = F)
    }
  }
  
  #generate output file name
  in.file <- gsub("\\..*", "", in.file)
  
  #convert gene id from ensembl to entrez gene id or viceversa and export the result
  if (gene_id=="ensembl"){
    cat("\n")
    cat("Convertion from ensembl ID to entrez ID ... \n")
    genelist <- biomaRt::getBM(attributes = c("ensembl_gene_id", "entrezgene_id", "external_gene_name"),
                               filters = c("ensembl_gene_id"),
                               values = genelist,
                               mart=(biomaRt::useMart("ensembl", biomart.species)), verbose = FALSE)
    colnames(genelist) <- c("ensembl_gene_id", "entrezgene", "external_gene_name")
    cat("DONE\n")
    #cat("\n")
    #cat("Checking for corresponding genes in KEGG database ... \n")
    cat("n.", paste(sum(length(which(!is.na(genelist$entrezgene))))), "out of", paste(sum(length(which(!is.na(genelist$ensembl_gene_id))))), "genes have a corresponding gene in KEGG database.", sep= " ")
    wd <- getwd()
    write.table(x = genelist, file = paste(in.file,"_converted.txt",sep=""), row.names = F)
    cat("\n")
    cat("Gene list exported!")
    return(genelist)
  }else{
    cat("\n")
    cat("Convertion from entrez ID to ensembl ID ... \n")
    genelist <- biomaRt::getBM(attributes = c("entrezgene_id", "ensembl_gene_id", "external_gene_name"),
                               filters = c("entrezgene_id"),
                               values = genelist,
                               mart=(biomaRt::useMart("ensembl", biomart.species)), verbose = FALSE)
    colnames(genelist) <- c("entrezgene", "ensembl_gene_id", "external_gene_name")
    cat("DONE\n")
    #cat("Checking for corresponding genes in KEGG database ... \n")
    cat("n.", paste(sum(length(which(!is.na(genelist$ensembl_gene_id))))), "out of", paste(sum(length(which(!is.na(genelist$entrezgene))))), "genes have a corresponding gene in KEGG database.", sep= " ")
    write.table(x = genelist, file = paste(in.file,"_converted.txt",sep=""), row.names = F)
    cat("\n")
    cat("Gene list exported!")
    return(genelist)
  }
}
