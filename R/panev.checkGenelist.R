#Script: panev.checkGenelist.R
#License: GPLv3 or later
#Modification date: 2019-05-13
#Written by: Valentino Palombo
#Contact: valentino.palombo@gmail.com
#Description: Check the genelist before run main PANEV script

panev.checkGenelist <- function(
  genelist=NULL, 
  type=NULL
)
{
  #check the arguments in the genelist
  if (is.null(genelist)){
    cat("Gene list not specified. Please check! \n")
    stop("Exit",call. = F)
  }else{
    cat("Gene list specified... ")
  }
  
  options(warn=-1)
  if (type == "gene"){
    if (colnames(genelist) %in% c("entrezgene", "ensembl_gene_id", "external_gene_name")){
      cat("and correct! \n")
    }else{
      cat("but is incorrect. Please check! \n")
      cat("Remember your gene list must have three columns with 'ensembl id', 'entrez id' and 'gene symbol', named 'ensembl_gene_id', 'entrezgene' and 'external_gene_name' respectively \n")
      stop("Exit",call. = F)
    }
  }
  if (type == "expr"){
    if ("ensembl_gene_id" %in% colnames(genelist) & "entrezgene" %in% colnames(genelist) &
        "FC" %in% colnames(genelist) & "pvalue" %in% colnames(genelist) & "external_gene_name" %in% colnames(genelist)){
      cat("and correct! \n")
    }else{
      cat("but is incorrect. Please check! \n")
      cat("Remember your gene list must have five columns with 'gene ensembl id', 'entrez gene id', 'gene symbol', 'Fold Change' and 'p-value', named 'ensembl_gene_id', 'entrezgene', 'external_gene_name', 'FC' and 'pvalue' respectively \n")
      stop("Exit",call. = F)
    }
  }
  options(warn=0)
}
