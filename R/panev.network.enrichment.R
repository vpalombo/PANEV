#' @title Enrichment analysis considering PANEV result as a background
#' @description The function perform an enrichment analysis considering as backgroud the union of pathways investigated by PANEV for generate network result.The result help to interpret the PANEV output allowing to identify the pathways most strongly associated with the input list of genes.
#' @usage panev.network.enrichment(in.file, out.file = "PANEV_enrich", species = NULL, FL = NULL, levels = 2)
#' @param in.file Name of input file (with extension) containing the gene list of interest. The file \bold{must} contain three columns labelled as \emph{'ensembl_gene_id'}, \emph{'entrezgene'} and \emph{'external_gene_name'}, respectively. The file \bold{must} rely in the working directory. The handy function \code{\link[PANEV]{panev.dataPreparation}} could be used to create a properly formatted input file from a single gene list.
#' @param out.file Name of the folder where the results will be stored and of the output diagram file (default =  'PANEV_enrich').
#' @param species The code of your species of interest. The correct code can get among the list of those available in KEGG with the handy function \code{\link[PANEV]{panev.speciesCode}}.
#' @param FL A list of pathways of first level to investigate. The list of all available pathways can get with the \code{\link[PANEV]{panev.pathList}} function. 
#' @param levels The number of levels of interactions (from \emph{1} to \emph{n}) investingated (default = 2).
#' @details This function is based on \code{\link[KEGGREST]{keggList}} and \code{\link[KEGGREST]{keggLink}} functions of \pkg{KEGGREST} package (\url{http://bioconductor.org/packages/release/bioc/html/KEGGREST.html}.
#' @details The enrichment analysis is based on \code{\link[bc3net]{enrichment}} function of \pkg{bc3net} package (\url{https://cran.r-project.org/web/packages/bc3net/index.html}.
#' @return A \emph{<out.file>_enrichment.txt} file containing the results of functional enrichment analysis, based on a one-sided Fisher's exact test (hypergeometric test) and considering the PANEV network as a background.
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
#' #Overview on input data 
#' panev.network.enrichment(in.file = in.file, 
#'               out.file = out.file, 
#'               species = species, 
#'               FL = FL, 
#'               levels = 3)
#' 

#Script: panev.network.enrichment.R
#License: GPLv3 or later
#Modification date: 2019-10-17
#Written by: Valentino Palombo
#Contact: valentino.palombo@gmail.com
#Description: Descriptive statistics of data (i.e. Enrichemnt analysis, Genes per pathways, Pathways for genes)

panev.network.enrichment <- function(
  in.file=NULL,
  out.file="PANEV_enrich",
  species=NULL,
  FL = NULL,
  levels = 2
)
{
  #import genelist 
  in.file.name <- paste(in.file)
  if (file.exists(in.file.name)){
    genelist <<- read.table(in.file.name, header=T, stringsAsFactors = F)
    cat("Input file imported! \n")
    cat("\n")
    #check the genelist arguments
    panev.checkGenelist(genelist = genelist, type="gene")
    #delete gene(s) with no entrez gene id
    genelist <- genelist[!is.na(genelist$entrezgene),]
    #add esembl to no symbol genes
    genelist <<- panev.addExternalGeneName(genelist = genelist)
  }else{
    cat("Input file not found. Please check! \n")
    stop("Exit",call. = F)
  }
  #check the path arguments
  panev.checkPath(species, FL)
  cat("\n")
  cat("\nPrerequisite check passed!\n\nPANEV is running ... \n")
  cat(" Please wait... It could be a while depending on the number of pathways and levels required! \n")
  #set the levels to investigate
  level <- 1:levels
  #download gene/path id table 
  panev.downloadGenePathTable(species = species)
  #add gene/path relations to genelist
  genelist$KEGGgenes <- paste(paste(species),":",genelist$entrezgene, sep="")
  genelist <- merge(x = genelist, y = KEGGgenes, by.x = "KEGGgenes", by.y = "V1", all=F)
  #download pathways description
  panev.downloadPathDescription(species = species)

  #prepare FL to query on KEGG
  query <- as.data.frame(gsub("map", species, FL), stringsAsFactors = F)
  colnames(query) <- "path_ID"
  #prepare the storage object
  nextlevel <- NULL
  if (length(level)>1){
    #path(s) excluded for further levels
    exclude <- query
    for (j in 2:length(level)){
      for (k in 1:length(query[,1])){
        out_nodes <- suppressMessages(data.frame(do.call(rbind, (XML::getNodeSet((XML::xmlParse(KEGGREST::keggGet(query[k,1], "kgml"))), "//entry[@type='map']", fun=XML::xmlToList)))))
        lks <- as.data.frame(t(as.data.frame(out_nodes$.attrs)), stringsAsFactors = F)
        lks <- unique(as.data.frame(paste(lks$name), stringsAsFactors = F))
        colnames(lks) <- "path_ID"
        lks <- as.data.frame(subset(lks, !(lks$path_ID %in% query$path_ID)), stringsAsFactors = F)
        nextlevel <- rbind(nextlevel, lks)
      } 
      nextlevel <- unique(nextlevel)
      #grab only path with species specified
      nextlevel <- nextlevel[grep(species, nextlevel$path_ID, ignore.case = TRUE), ]
      #grab path different from previous level
      nextlevel <- nextlevel[(!(nextlevel %in% exclude$path_ID))]
      #path excluded in the next levels analysis
      exclude <- rbind(exclude, data.frame(path_ID = nextlevel))
      query <- data.frame(path_ID = nextlevel, stringsAsFactors = F) 
    }
  exclude <- unique(exclude)
  exclude <- exclude$path_ID
  cat("\nEnrichment analysis started ... ")
  #perform the enrichment analysis
  genes <- as.vector(unique(genelist$KEGGgenes))
  KEGGgenes2 <- as.data.frame(KEGGgenes)
  KEGGgenes2$V1 <- as.character(KEGGgenes2$V1)
  pathquery <- gsub("map", species, exclude)
  KEGGgenes2 <- KEGGgenes2[KEGGgenes2$rn %in% exclude,]
  reference <- unique(paste(KEGGgenes2$V1))
  pathtofind <- unique(exclude)
  genesets <- list()
  for (i in 1:length(pathtofind)) {
    genesets[[paste(pathtofind[i])]] <- subset(KEGGgenes$V1, KEGGgenes$rn==pathtofind[i])
  }
  enrich <- bc3net::enrichment(genes, reference, genesets, adj = "fdr", verbose = FALSE)
  enrich <- merge(x = enrich, y = KEGGpath, by.x = "TermID", by.y = "rn", all=F)
  colnames(enrich) <- c("pathway_ID", "n_genes", "all_genes", "pvalue", "padj", "pathway_name")
  enrich <- enrich[order(enrich$padj, enrich$pvalue),]
  write.table(x = enrich, file = paste(out.file,"_enrichment.txt",sep=""), row.names = F)
  cat("\n")
  cat("Enrichment results exported! \n")
  }else{
    exclude <- FL
    exclude <- gsub("map", species, exclude)
    cat("\nEnrichment analysis started ... ")
    #perform the enrichment analysis
    genes <- as.vector(unique(genelist$KEGGgenes))
    KEGGgenes2 <- as.data.frame(KEGGgenes)
    KEGGgenes2$V1 <- as.character(KEGGgenes2$V1)
    pathquery <- gsub("map", species, exclude)
    KEGGgenes2 <- KEGGgenes2[KEGGgenes2$rn %in% exclude,]
    reference <- unique(paste(KEGGgenes2$V1))
    pathtofind <- unique(exclude)
    genesets <- list()
    for (i in 1:length(pathtofind)) {
      genesets[[paste(pathtofind[i])]] <- subset(KEGGgenes$V1, KEGGgenes$rn==pathtofind[i])
    }
    enrich <- bc3net::enrichment(genes, reference, genesets, adj = "fdr", verbose = FALSE)
    enrich <- merge(x = enrich, y = KEGGpath, by.x = "TermID", by.y = "rn", all=F)
    colnames(enrich) <- c("pathway_ID", "n_genes", "all_genes", "pvalue", "padj", "pathway_name")
    enrich <- enrich[order(enrich$padj, enrich$pvalue),]
    write.table(x = enrich, file = paste(out.file,"_enrichment.txt",sep=""), row.names = F)
    cat("\n")
    cat("and results exported! \n")
  }
}
  