#' @title Data overview and enrichment analysis
#' @description The function takes as input the principal descriptive information on given gene list. In particular: number of gene(s) per single pathway, number of pathway(s) per single gene and enrichment analysis. When biological assumptions are broad or not strictly defined, pathways with more occurrences may be principal candidate pathways (FL) to be investigated with PANEV. The results obtained are based on KEGG Pathway database information.
#' @usage panev.stats.enrichment(in.file, out.file = "PANEV_enrich", species = NULL)
#' @param in.file Name of input file (with extension). The file \bold{must} contain three columns labelled \emph{'ensembl_gene_id'}, \emph{'entrezgene'} and \emph{'external_gene_name'} respectively. The file \bold{must} rely in the working directory. The handy function \code{\link[PANEV]{panev.dataPreparation}} could be used to create a properly formatted input file from a gene list.
#' @param out.file The specific name of generated files (default = 'PANEV_enrich').
#' @param species The species code of interest. You can get the correct code among the list of those available in KEGG with the \code{\link[PANEV]{panev.speciesCode}} function.
#' @details This function is based on \code{\link[KEGGREST]{keggList}} and \code{\link[KEGGREST]{keggLink}} functions of \pkg{KEGGREST} package (\url{http://bioconductor.org/packages/release/bioc/html/KEGGREST.html}.
#' @details The enrichment analysis is based on \code{\link[bc3net]{enrichment}} function of \pkg{bc3net} package (\url{https://cran.r-project.org/web/packages/bc3net/index.html}.
#' @return A \emph{<out.file>_enrichment.txt} file containing the results of functional enrichment analysis, based on a one-sided Fisher's exact test (hypergeometric test).
#' @return A \emph{<out.file>_GxP.txt} file containing the numbers of gene(s) detected for single pathway.
#' @return A \emph{<out.file>_PxG.txt} file containing the numbers of pathway(s) detected for single gene.
#' @author Valentino Palombo (\email{valentino.palombo@gmail.com})
#' @references Tenenbaum D (2017). KEGGREST: Client-side REST access to KEGG. R package version 1.16.1. 
#' @references Simoes R de M, Emmert-Streib F (2012). Bagging statistical network inference from large-scale gene expression data. PLOS ONE; 7: e33624. doi:10.1371/journal.pone.0033624
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
#' #Look for the specie code matching the search string 
#' list <- panev.speciesCode(string = "bos")
#' species <- as.character(list[1,2]) # bta
#' 
#' #Run the PANEV enrichment function  
#' panev.stats.enrichment(in.file = "ensembl_genelist_converted.txt", 
#'                      out.file = "example", 
#'                      species= species)
#' 

#Script: panev.stats.enrichment.R
#License: GPLv3 or later
#Modification date: 2019-05-13
#Written by: Valentino Palombo
#Contact: valentino.palombo@gmail.com
#Description: Descriptive statistics of data (i.e. Enrichemnt analysis, Genes per pathways, Pathways for genes)

panev.stats.enrichment <- function(
  in.file=NULL,
  out.file="PANEV_enrich",
  species=NULL
  )
{
  #cat("Arguments checking ... \n")
  #import genelist 
  in.file.name <- paste(in.file)
  if (file.exists(in.file.name)){
    genelist <<- read.table(in.file.name, header=T, stringsAsFactors = F)
    cat("Input file is imported! \n")
    cat("\n")
    #add esembl to no symbol genes
    genelist <<- panev.addExternalGeneName(genelist = genelist)
  }else{
    cat("Input file not found. Please check! \n")
    stop("Exit",call. = F)
  }
  #check genelist 
  panev.checkGenelist(genelist = genelist, type="gene")
  #check specie code
  panev.checkSpecies(species)
  
  # Function #
  cat("\nEnrichment analysis started ... ")
  #download genes list available with the relative root pathway(s)
  KEGGgenes <- data.frame(cbind(rn=names(KEGGREST::keggLink(paste(species),"pathway")), 
                                V1=KEGGREST::keggLink(paste(species),"pathway")), stringsAsFactors = F, row.names = NULL)
  genelist$KEGGgenes <- paste(paste(species),":",genelist$entrezgene, sep="")
  genelist <- merge(x = genelist, y = KEGGgenes, by.x = "KEGGgenes", by.y = "V1", all=F)
  #download the available KEGG pathways list  
  KEGGpath <- data.table::setDT(as.data.frame(cbind(KEGGREST::keggList("pathway", species))), keep.rownames = T)
  orgs <- as.data.frame(KEGGREST::keggList("organism"))[2:3]
  string <- paste(" -", paste(orgs[which(orgs$organism == species),2]))    
  KEGGpath$V1 <- stringr::str_replace_all(KEGGpath$V1, stringr::fixed(string), "")
  #perform and enrichment analysis
  genes <- as.vector(unique(genelist$KEGGgenes))
  KEGGgenes2 <- as.data.frame(KEGGgenes)
  KEGGgenes2$V1 <- as.character(KEGGgenes2$V1)
  reference <- unique(paste(KEGGgenes2$V1))
  pathtofind <- unique(KEGGgenes2$rn)
  genesets <- list()
  for (i in 1:length(pathtofind)) {
    genesets[[paste(pathtofind[i])]] <- subset(KEGGgenes2$V1, KEGGgenes2$rn==pathtofind[i])
  }
  enrich <- bc3net::enrichment(genes, reference, genesets, adj = "fdr", verbose = FALSE)
  enrich <- merge(x = enrich, y = KEGGpath, by.x = "TermID", by.y = "rn", all=F)
  colnames(enrich) <- c("pathway_ID", "n_genes", "all_genes", "pvalue", "padj", "pathway_name")
  enrich <- enrich[order(enrich$padj, enrich$pvalue),]
  write.table(x = enrich, file = paste(out.file,"_enrichment.txt",sep=""), row.names = F)
  cat("\n")
  cat("and results exported! \n")
  #create and save the gene(s) per pathway(s) frequency table 
  freqGxP <- merge(x = as.data.frame(table(genelist$rn)), y = KEGGpath, by.x = "Var1", by.y = "rn", all=F)
  freqGxP <- data.frame(freqGxP[order(freqGxP$Freq, decreasing = T), ], row.names = NULL)
  freqGxP <- freqGxP[,c(2,3,1)]
  colnames(freqGxP) <- c("n_genes", "pathway_name", "pathway_ID")
  write.table(x = freqGxP, file = paste(out.file,"_GxP.txt",sep=""), row.names = F)
  cat("\n")
  cat("Gene per pathway(s) table created and exported! \n")
  #create and save the pathway(s) per gene(s) frequency table
  freqPxG <- merge(x = data.frame(table(genelist$entrezgene)), y = genelist[c("entrezgene","ensembl_gene_id", "external_gene_name")], by.x = "Var1", by.y = "entrezgene", all=F)
  freqPxG <- unique(data.frame(freqPxG[order(freqPxG$Freq, decreasing = T), ], row.names = NULL))
  freqPxG <- freqPxG[,c(2,1,3,4)]
  colnames(freqPxG) <- c("n_pathways", "entrez_gene_id", "ensembl_gene_id","gene_symbol")
  write.table(x = freqPxG, file = paste(out.file,"_PxG.txt",sep=""), row.names = F)
  cat("\n")
  cat("Pathway per gene(s) table created and exported! \n")
}