#' @title PANEV on gene expression datasets
#' @description The function visualizes PANEV results from expression gene datasets, considering any interactions among selected pathways. For each pathway the diagram helps to visualize associated genes that overcome a user defined \emph{p-value} threshold. The color is based on the fold change (\emph{FC}) value.
#' @usage panev.exprnetwork(in.file = NULL, path.file = NULL, out.file = "PANEV_expr", 
#'                       species = NULL, pvalue = 0.05)
#' @param in.file Name of input file (with extension). The input file containing the differentially expressed genes (DEG) analysis results. The file \bold{must} contain five columns labeled \emph{'ensembl_gene_id'}, \emph{'entrezgene'}, \emph{'external_gene_name'}, \emph{'FC'} and \emph{'pvalue'}. The \code{\link[PANEV]{panev.exprdataPreparation}} function can be used to create the properly formatted input file from a gene expression list. The input file \bold{must} rely in working directory (default = NULL).
#' @param path.file Name of input file (with extension). The input file containing two columns labeled \emph{path_ID} and \emph{value}, containing respectively the ID of pathways of interest (\code{\link[PANEV]{panev.pathList}} function could be used to check if they are available), and relative biological estimated score value (e.g. flux value obtained with Dinamic Impact Approach analysis developed and described by Bionaz et al., 2012). The input file \bold{must} rely in working directory (default = NULL).
#' @param out.file Name of output diagram file and folder where results will be stored (default = "PANEV_expr").
#' @param species The KEGG organism code. The correct code for the species of interest can be retrieved with the handy \code{\link[PANEV]{panev.speciesCode}} function (default = NULL).
#' @param pvalue The p-value cut-off for Differentially Expressed (DE) analysis (default = \emph{0.05}).
#' @details This function is based on the main KEGGREST query function \code{\link[KEGGREST]{keggLink}}.
#' @return A \emph{<out.file>.html} file with the diagram visualization of PANEV results is created and stored in a folder named as \emph{PANEV_exprRESULTS_<out.name>}, created in the work directory.
#' @author Valentino Palombo (\email{valentino.palombo@gmail.com})
#' @references Bionaz, M., K. Periasamy, S.L. Rodriguez-Zas, W.L. Hurley, and J.J. Loor. 2012. A novel dynamic impact approach (DIA) for functional analysis of time-course omics studies: validation using the bovine mammary transcriptome. PLos One 7:e32455. doi:10.1371/journal.pone.0032455. 
#' @references Tenenbaum D (2017). KEGGREST: Client-side REST access to KEGG. R package version 1.16.1. 
#' @references Thieurmel B (2016). visNetwork: Network Visualization using 'vis.js' Library. R package version 2.0.3. https://CRAN.R-project.org/package=visNetwork
#' @examples ##### EXAMPLES CODE #####
#' #Copy the example files in the current working directory
#' panev.example()
#' 
#' #Parameters
#' in.file = "exprdata.txt"
#' path.file  = "expr_listPath.txt"
#' out.file = "expression_data"
#' pvalue = 0.05
#' 
#' #Look for the specie code matching the search string 
#' list <- panev.speciesCode(string = "scrofa")
#' species = as.character(list[1,2]) # ssc
#' 
#' # Run the PANEV function
#' panev.exprnetwork(in.file = in.file, 
#'                   path.file = path.file, 
#'                   out.file = out.file, 
#'                   species = species, 
#'                   pvalue = pvalue)
#' 

#Script: panev.exprnetwork.R
#License: GPLv3 or later
#Modification date: 2019-10-15
#Written by: Valentino Palombo
#Contact: valentino.palombo@gmail.com
#Description: PANEV diagram visualization of expression datasets.

panev.exprnetwork <- function(
  in.file = NULL,
  path.file = NULL,
  out.file = "PANEV_expr",
  species = NULL,
  pvalue = 0.05
)
{
  #import genelist 
  in.file.name <- paste(in.file)
  if (file.exists(in.file.name)){
    genelist <<- read.table(in.file.name, header=T, stringsAsFactors = F)
    cat("Input file imported! \n")
    #add esembl to no symbol genes
    genelist <<- panev.addExternalGeneName(genelist = genelist)
  }else{
    cat("Input file not found. Please check! \n")
    stop("Exit",call. = F)
  }
  #import path info
  in.file.name <- paste(path.file)
  if (file.exists(in.file.name)){
    listPath <<- read.table(in.file.name, header=T, stringsAsFactors = F)
    cat("\n")
    cat("Pathway input file imported! \n")
  }else{
    cat("Pathway input file not found. Please check! \n")
    stop("Exit",call. = F)
  }
  cat("\n")
  #Check the genes
  panev.checkGenelist(genelist = genelist, type="expr")
  #Check the path arguments
  if ("path_ID" %in% colnames(listPath) & "value" %in% colnames(listPath)){
    cat("\nYour path list colnames are correct! \n")
  }else{
    cat("\nYour path list colnames are incorrect. Please check! \n")
    cat("Remember your path list must have two columns with a list of pathways of interest and an expression value related, named 'path_ID' and 'value' respectively \n")
    stop("Exit",call. = F)
  }
  #read the list of pathway(s) to investigate
  FL <<- listPath$path_ID
  panev.checkPath(FL = FL, species = species)
  cat("\n")
  cat("\nPrerequisite check passed!\n\nPANEV is running ... \n")
  cat(" Please wait... It could take a while depending on the number of pathways required! \n")
  #pvalue cut-off
  ngenelist <- nrow(genelist)
  genelist <- genelist[genelist$pvalue < pvalue,]
  cat("n.", nrow(genelist), "of", ngenelist, "genes passed the p-value filtering.\n", sep= " ")
  #download pathways id table (with related genes) and add to genelist
  KEGGgenes <- panev.downloadGenePathTable(species = species)
  #add gene/path edges to genelist
  genelist$KEGGgenes <- paste(paste(species),":",genelist$entrezgene, sep="")
  genelist <- merge(x = genelist, y = KEGGgenes, by.x = "KEGGgenes", by.y = "V1", all=F)
  #download pathways description  
  panev.downloadPathDescription(species = species)
  #change the name of columns
  colnames(KEGGpath) <- c("rn", "path")
  #select di FL pathway(s) provided
  FL <<- data.frame(gsub("map", species, listPath$path_ID), stringsAsFactors = F)
  colnames(FL) <<- "path_ID" # check
  #extract the list of pathways to analyse
  list_pathway <<- as.data.frame(table(genelist$rn), stringsAsFactors = F)
  list_pathway <- list_pathway[list_pathway$Var1 %in% FL$path_ID,]    
  if(length(list_pathway$Var1)==0){
    cat("Your pathway(s) of interest is/are not present among the list of pathways available for analysis. \n")
    cat("This may depend on your pvalue cut-off (default =  0.05), \n")
    cat("or simply because your pathway(s) of interest is/are not represented by your DEGs list")
    stop("Exit",call. = F)
  }else{
    #prepare FL to query on KEGG
    query <- as.data.frame(list_pathway$Var1, stringsAsFactors = F)
    colnames(query) <- "path_ID"
    nextlevel <- NULL
    #query on KEGG
    for (k in 1:length(query[,1])){
      out_nodes <- suppressMessages(data.frame(do.call(rbind, (XML::getNodeSet((XML::xmlParse(KEGGREST::keggGet(query[k,1], "kgml"))), "//entry[@type='map']", fun=XML::xmlToList)))))
      lks <- as.data.frame(t(as.data.frame(out_nodes$.attrs)), stringsAsFactors = F)
      lks <- unique(as.data.frame(paste(lks$name), stringsAsFactors = F))
      colnames(lks) <- "path_ID"
      lks <- as.data.frame(subset(lks, lks$path_ID %in% paste(query$path_ID)), stringsAsFactors = F)
      lks$to <- lks$path_ID
      lks$from <- paste(query[k,1])
      lks <- lks[grep(species, lks$to), ]
      nextlevel <- rbind(nextlevel, lks)
    }
    #select the interactions and store the results
    int <- as.data.frame(cbind(nextlevel$from, nextlevel$to), stringsAsFactors = F, row.names = F)
    colnames(int) <- c("from", "to")
    int <- as.data.frame(dplyr::filter(int, int$from %in% query$path_ID), stringsAsFactors = F)
    int <- as.data.frame(dplyr::filter(int, int$to %in% query$path_ID), stringsAsFactors = F)
    int <- int[nextlevel$from != nextlevel$to,]
    results <- as.data.frame(dplyr::filter(listPath, gsub("map", species, listPath$path_ID) %in% query$path_ID), stringsAsFactors = F)
    results$path_ID<- gsub("map", species, results$path_ID)
    colnames(results) <- c("V1", "value")
    results$shape <- "diamond"
    FCdata <- genelist[genelist$rn %in% unique(nextlevel$from),]
    int2 <- as.data.frame(cbind(FCdata$external_gene_name, FCdata$rn), stringsAsFactors = F)
    colnames(int2) <- c("from", "to")
    int <- as.data.frame(rbind(int2, int))
    int <- unique(int)
    FCdata <- data.frame(cbind(paste(FCdata$external_gene_name), as.numeric(FCdata$FC)))
    colnames(FCdata) <- c("V1", "value")
    FCdata$V1 <- paste(FCdata$V1)
    FCdata$value <- as.numeric(paste(FCdata$value))
    FCdata$shape <- "circle"
    FCdata <- unique(FCdata)
    results <- as.data.frame(rbind(FCdata, results))
    colnames(results) <- c("name", "value", "shape")
    label <- results
    label$name <- paste(label$name)
    #merge with path description
    label$name[label$name %in% KEGGpath$rn] <- KEGGpath$path[KEGGpath$rn %in% label$name]
    label$value <- as.numeric(paste(label$value))
    label$legend <- c(1:length(label$name))
    label$value[label$shape == "diamond"] <- label$value[label$shape == "diamond"]/max(abs(label$value[label$shape == "diamond"]))
    label$value[label$shape == "circle"] <- label$value[label$shape == "circle"]/max(abs(label$value[label$shape == "circle"]))
    label$type <- "1"
    label$type <- ifelse(label$shape=="circle"& label$value>=0.0 & label$value<0.25,"low upregulated gene", label$type)
    label$type <- ifelse(label$shape=="circle"& label$value>=0.25 & label$value<0.50,"moderate upregulated gene", label$type)
    label$type <- ifelse(label$shape=="circle"& label$value>=0.50 & label$value<0.75,"high upregulated gene", label$type)
    label$type <- ifelse(label$shape=="circle"& label$value>=0.75 & label$value<=1,"strong upregulated gene", label$type)
    label$type <- ifelse(label$shape=="circle"& label$value<=0.0 & label$value>(-0.25),"low downregulated gene", label$type)
    label$type <- ifelse(label$shape=="circle"& label$value<=(-0.25) & label$value>(-0.50),"moderate downregulated gene", label$type)
    label$type <- ifelse(label$shape=="circle"& label$value<=(-0.50) & label$value>(-0.75),"high downregulated gene", label$type)
    label$type <- ifelse(label$shape=="circle"& label$value<=(-0.75) & label$value>=(-1),"strong downregulated gene", label$type)
    label$type <- ifelse(label$shape=="diamond"& label$value>=0.0 & label$value<0.25,"low upregulated pathway", label$type)
    label$type <- ifelse(label$shape=="diamond"& label$value>=0.25 & label$value<0.50,"moderate upregulated pathway", label$type)
    label$type <- ifelse(label$shape=="diamond"& label$value>=0.50 & label$value<0.75,"high upregulated pathway", label$type)
    label$type <- ifelse(label$shape=="diamond"& label$value>=0.75 & label$value<=1,"strong upregulated pathway", label$type)
    label$type <- ifelse(label$shape=="diamond"& label$value<=0.0 & label$value>(-0.25),"low downregulated pathway", label$type)
    label$type <- ifelse(label$shape=="diamond"& label$value<=(-0.25) & label$value>(-0.50),"moderate downregulated pathway", label$type)
    label$type <- ifelse(label$shape=="diamond"& label$value<=(-0.50) & label$value>(-0.75),"high downregulated pathway", label$type)
    label$type <- ifelse(label$shape=="diamond"& label$value<=(-0.75) & label$value>=(-1),"strong downregulated pathway", label$type)
    #create nodes
    ndf <- data.frame(
      id = label$legend,
      label = paste(label$name),
      shape = label$shape,
      group = label$type)
    #create edges 
    dm <- int
    #merge with path description
    for (i in 1:length(KEGGpath$rn)){
      dm$from[paste(dm$from) == paste(KEGGpath$rn[i])] <- paste(KEGGpath$path[i])
      dm$to[paste(dm$to) == paste(KEGGpath$rn[i])] <- paste(KEGGpath$path[i])
    }
    #add id number
    for (i in 1:length(label$legend)){
      dm$from[paste(dm$from) == paste(label$name[i])] <- paste(label$legend[i])
      dm$to[paste(dm$to) == paste(label$name[i])] <- paste(label$legend[i])
    }
    #create the diagram edges
    edf <- data.frame(
      from = dm$from,
      to = dm$to)
    #legend color
    red <- colorRampPalette(RColorBrewer::brewer.pal(4,"Reds")) (4)
    green <- colorRampPalette(RColorBrewer::brewer.pal(4,"Greens")) (4)
    #create the diagram and export a html file
    graph <- visNetwork::visNetwork(ndf, edf, main = paste("PANEV visualization", out.file, sep=" "))  
    graph <- visNetwork::visGroups(graph, groupname = "low upregulated gene", color = paste(red[1]), shape="ellipse")  
    graph <- visNetwork::visGroups(graph, groupname = "moderate upregulated gene", color = paste(red[2]), shape="ellipse") 
    graph <- visNetwork::visGroups(graph, groupname = "high upregulated gene", color = paste(red[3]), shape="ellipse") 
    graph <- visNetwork::visGroups(graph, groupname = "strong upregulated gene", color = paste(red[4]), shape="ellipse") 
    graph <- visNetwork::visGroups(graph, groupname = "low downregulated gene", color = paste(green[1]), shape="ellipse") 
    graph <- visNetwork::visGroups(graph, groupname = "moderate downregulated gene", color = paste(green[2]), shape="ellipse") 
    graph <- visNetwork::visGroups(graph, groupname = "high downregulated gene", color = paste(green[3]), shape="ellipse") 
    graph <- visNetwork::visGroups(graph, groupname = "strong downregulated gene", color = paste(green[4]), shape="ellipse") 
    graph <- visNetwork::visGroups(graph, groupname = "low upregulated pathway", color = paste(red[1]), shape="diamond")  
    graph <- visNetwork::visGroups(graph, groupname = "moderate upregulated pathway", color = paste(red[2]), shape="diamond")  
    graph <- visNetwork::visGroups(graph, groupname = "high upregulated pathway", color = paste(red[3]), shape="diamond")  
    graph <- visNetwork::visGroups(graph, groupname = "strong upregulated pathway", color = paste(red[4]), shape="diamond")  
    graph <- visNetwork::visGroups(graph, groupname = "low downregulated pathway", color = paste(green[1]), shape="diamond")  
    graph <- visNetwork::visGroups(graph, groupname = "moderate downregulated pathway", color = paste(green[2]), shape="diamond")  
    graph <- visNetwork::visGroups(graph, groupname = "high downregulated pathway", color = paste(green[3]), shape="diamond")  
    graph <- visNetwork::visGroups(graph, groupname = "strong downregulated pathway", color = paste(green[4]), shape="diamond")  
    graph <- visNetwork::visOptions(graph, highlightNearest = TRUE, nodesIdSelection = TRUE)
    graph <- visNetwork::visLegend(graph, useGroups = T)
    #export the results
    RES <- paste("PANEV_exprRESULTS_",out.file,sep="")
    if (file.exists(RES)) {
      RES <- paste(RES,"/",sep="")
      wd <- getwd()
      visNetwork::visSave(graph, file = paste(wd, "/", RES,out.file,".html", sep=""))
    }else{
      dir.create(RES)
      RES <- paste(RES,"/",sep="")
      wd <- getwd()
      visNetwork::visSave(graph, file = paste(wd, "/",RES,out.file,".html", sep=""))
    }
    #clear the workspace
    panev.cleanEnvir(KEGGpath)
    panev.cleanEnvir(KEGGgenes)
    panev.cleanEnvir(list_pathway)
    panev.cleanEnvir(FL)
    cat("\n")
    cat("Well done! Diagram visualization was created and exported. \n")
  }
}  
