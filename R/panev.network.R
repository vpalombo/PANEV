#' @title Pathways Network Visualization (PANEV)
#' @description The function performs a pathway analysis taking into account both upstream and downstream dependent network of functional related genes, from 2 to \emph{n} degrees of interaction. This function generates \emph{n} '\emph{.txt}' files, one for each level of interaction, containing the candidate genes and the related pathways highlighted. Along with the tabular format results, the function gets also the diagram visualization of PANEV results, saved in an interactive '\emph{.html}' file.
#' @usage panev.network(in.file = NULL, out.file = "PANEV_gene", species = NULL, 
#'                   FL = NULL, levels = 2)
#' @param in.file Name of input file (with extension) containing the gene list of interest. The file \bold{must} contain three columns labelled as \emph{'ensembl_gene_id'}, \emph{'entrezgene'} and \emph{'external_gene_name'}, respectively. The file \bold{must} rely in the working directory. The handy function \code{\link[PANEV]{panev.dataPreparation}} could be used to create a properly formatted input file from a single gene list.
#' @param out.file Name of the folder where the results will be stored and of the output diagram file (default =  'PANEV_gene').
#' @param species The code of your species of interest. The correct code can get among the list of those available in KEGG with the handy function \code{\link[PANEV]{panev.speciesCode}}.
#' @param FL A list of pathways of first level to investigate. The list of all available pathways can get with the \code{\link[PANEV]{panev.pathList}} function. 
#' @param levels The number of levels of interactions (from \emph{1} to \emph{n}) investingated (default = 2).
#' @details For each level required for the investigation, a single \emph{.txt} file will be generated, containing the genes and the related pathways for a specific level of interaction. The function get also the diagram visualization of results. This function is based on the main \pkg{KEGGREST} package (\url{https://bioconductor.org/packages/release/bioc/html/KEGGREST.html}).
#' @return 
#' \itemize{
#' \item{The script generates from \emph{1} to \emph{n} '\emph{txt}' files, (with \emph{n} equal to the levels of interaction required for the investigation). 
#' Each file is named \emph{(n)}Lgenes, based on \emph{n} level analysed 
#' and contains the genes falling inside the pathways of the specific degree of interaction investigated.
#' Each file contains five columns with \emph{ensembl gene}, \emph{entrez gene}, \emph{gene name}, \emph{path description} and \emph{path id}. 
#' The files are stored in a folder named as \emph{PANEV_RESULTS_<out.name>}, created in the work directory.}
#' \item{A \emph{<out.file>.html} file with the diagram visualization of PANEV results. The file is stored in the same folder of \emph{.txt} files}
#' }
#' @author Valentino Palombo (\email{valentino.palombo@gmail.com})
#' @references Tenenbaum D (2017). KEGGREST: Client-side REST access to KEGG. R package version 1.16.1. 
#' @references Thieurmel B (2016). visNetwork: Network Visualization using 'vis.js' Library. R package version 2.0.3.  https://CRAN.R-project.org/package=visNetwork
#' @examples ##### EXAMPLES CODE #####
#' #Copy the example files in the current working directory
#' panev.example()
#' 
#' #Parameters
#' in.file="data.txt"
#' out.file="example"
#' FL = c("path:map00061", "path:map00062", "path:map00071", "path:map00072")
#' levels = 2
#' 
#' #Look for the specie code matching the search string 
#' list <- panev.speciesCode(string = "bos")
#' species=as.character(list[1,2]) # bta 
#' 
#' # Run the PANEV function
#' panev.network(in.file = in.file, 
#'               out.file = out.file, 
#'               species = species, 
#'               FL = FL, 
#'               levels = levels)
#'            

#Script: panev.network.R
#License: GPLv3 or later
#Modification date: 2019-10-15
#Written by: Valentino Palombo
#Contact: valentino.palombo@gmail.com
#Description: Run PANEV and get diagram visualization.

panev.network <- function(
  in.file=NULL,
  out.file="PANEV_gene",
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
  #create a folder for storing results
  RES <- paste("PANEV_RESULTS_",out.file,sep="")
  if (file.exists(RES)) {
    RES <- paste(RES,"/",sep="")
    wd <- getwd()
  }else{
    dir.create(RES)
    RES <- paste(RES,"/",sep="")
    wd <- getwd()
  }
  #prepare FL to query on KEGG
  query <- as.data.frame(gsub("map", species, FL), stringsAsFactors = F)
  colnames(query) <- "path_ID"
  #prepare the storage object
  nextlevel <- NULL
  glFL <- as.data.frame(dplyr::filter(genelist, genelist$rn %in% query$path_ID), stringsAsFactors = F)
  FLgenes <- as.data.frame(glFL, stringsAsFactors = F)
  FLgenes$external_gene_name <- as.character(FLgenes$external_gene_name)
  FLgenes$ensembl_gene_id <- as.character(FLgenes$ensembl_gene_id) 
  FLgenes <- data.frame(cbind(ensembl_gene_id = paste(FLgenes$ensembl_gene_id), entrezgene = paste(FLgenes$entrezgene), external_gene_name = paste(FLgenes$external_gene_name), rn = paste(FLgenes$rn)), stringsAsFactors = F)
  FLgenes <- merge(FLgenes, KEGGpath, by.x = "rn", by.y = "rn")
  FLgenes <- data.frame(cbind(paste(FLgenes$ensembl_gene_id), paste(FLgenes$entrezgene), paste(FLgenes$external_gene_name), paste(FLgenes$V1), paste(FLgenes$rn)), stringsAsFactors = F)
  colnames(FLgenes) <- c("ensemblgene","entrezgene", "gene_name", "path_description","path_ID")
  exp <- paste(RES,"/1Lgenes.txt",sep="")
  write.table(FLgenes, exp, sep="\t", row.names = F)
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
      glNDI <- as.data.frame(dplyr::filter(genelist, genelist$rn %in% nextlevel), stringsAsFactors = F)
      glNDI <- merge(glNDI, KEGGpath, by.x = "rn", by.y = "rn")
      glNDI <- as.data.frame(cbind(paste(glNDI$ensembl_gene_id), glNDI$entrezgene, paste(glNDI$external_gene_name), glNDI$V1, glNDI$rn))
      colnames(glNDI) <- c("ensemblgene","entrezgene", "gene_name", "path_description","path_ID")
      names <- paste(j,"Lgenes.txt", sep="")
      exp <- paste(RES, names, sep="")
      write.table(glNDI, exp, sep="\t", row.names = F) 
      query <- data.frame(path_ID = nextlevel, stringsAsFactors = F) 
    }
  }
  cat("\n")
  cat("PANEV analysis completed and relative '.txt' files exported! \n")
  #set wd where PANEV outputs were stored
  setwd(RES)
  #upload the output file(s) generated before 
  panev.uploadPANEVscriptoutput(levels=levels)
  cat("\n")
  cat("Preparing PANEV diagram visualization! \n")
  cat(" Please wait... It could be a while depending on the number of pathways and levels required! \n")
  #find the interactions
  #create the species-specific pathway(s) list for querying the KEGG database
  if (length(level)==1){
    query <- as.data.frame(gsub("map", species, FL), stringsAsFactors = F)
    colnames(query)="path_ID"
    interactions <- NULL
    for (k in 1:length(query[,1])){
      out_nodes <- suppressMessages(data.frame(do.call(rbind, (XML::getNodeSet((XML::xmlParse(KEGGREST::keggGet(query[k,1], "kgml"))), "//entry[@type='map']", fun=XML::xmlToList)))))
      lks <- as.data.frame(t(as.data.frame(out_nodes$.attrs)), stringsAsFactors = F)
      lks <- unique(as.data.frame(paste(lks$name), stringsAsFactors = F))
      colnames(lks) <- "path_ID"
      lks <- as.data.frame(subset(lks, (lks$path_ID %in% query$path_ID)), stringsAsFactors = F)
      lks$to <- lks$path_ID
      lks$from <- query[k,1]
      lks <- lks[grep(species, lks$to), ]
      interactions <- rbind(interactions, lks)
    }
    interactions$path_ID <-NULL
    interactions <- unique(interactions)
    interactions <- interactions[interactions$from != interactions$to,]
    dm <- NULL
    dm <- rbind(dm, get(paste("the","1Lgenes",sep="")))
    dm <- as.data.frame(cbind(dm$gene_name, dm$path_ID), row.names = F, stringsAsFactors = F)
    colnames(dm) <- c("from", "to")
    int <- as.data.frame(cbind(interactions$from, interactions$to), row.names = F)
    colnames(int) <- c("from", "to")
    dm <- as.data.frame(rbind(dm, int))
    row.names(dm) <- NULL
    #extract the interactions uncovered to generate the diagram label
    label <- data.frame(name = unique(c(paste(dm$from), paste(dm$to))))
    label$num <- c(1:length(label$name))
    #create the attributes of diagram objects
    label$type <- 1
    add <- get(paste("the","1Lgenes",sep=""))
    label$type <- ifelse(label$name %in% add$gene_name, paste("1Lgenes", sep=""), label$type)
    label$type <- ifelse(label$name %in% add$path_ID, paste("1Lpath", sep=""), label$type)
    label$shape <- 1
    label$shape <- ifelse(label$name %in% add$gene_name, "circle", label$shape)
    label$shape <- ifelse(label$shape == "1", "diamond", label$shape)
    label$color <- "black"
    label$type <- ifelse(label$type == "1", "path without results", label$type)
    label$type <- ifelse(label$shape == "circle", "gene", label$type)
    label <- label[order(label$type, decreasing = F),]
    label <- label[order(label$shape, decreasing = F),]
    #create the diagram legend
    legenda <- data.frame(cbind(num = label$num, name = paste(label$name)), stringsAsFactors = F)
    legenda <- merge(legenda, KEGGpath, by.x = "name", by.y = "rn", all=T)
    legenda <- legenda[order(as.numeric(legenda$num)),]
    legenda <- subset(legenda, (!is.na(legenda$num)))
    legenda$name <- ifelse(is.na(legenda$V1), legenda$name, legenda$V1)
    #create diagram nodes
    ndf <- data.frame(
      id = label$num,
      label = paste(legenda$name),
      shape = label$shape,
      group = label$type)
    #extract the interactions uncovered 
    for (i in 1:length(label$num)){
      dm$from[paste(dm$from) == paste(label$name[i])] <- paste(label$num[i])
      dm$to[paste(dm$to) == paste(label$name[i])] <- paste(label$num[i])
    }
    #create the diagram edges
    edf <- data.frame(
      from = dm$from,
      to = dm$to)
    #legend color
    qual_col_pals = RColorBrewer::brewer.pal.info[RColorBrewer::brewer.pal.info$category == 'qual',]
    col_vector = unlist(mapply(RColorBrewer::brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
    #create the diagram and export a html file
    graph <- visNetwork::visNetwork(ndf, edf, main = paste("PANEV visualization", out.file, sep=" "))  
    graph <- visNetwork::visGroups(graph, groupname = "gene", color = paste(col_vector[1]), shape="dot")  
    graph <- visNetwork::visGroups(graph, groupname = "path without results", color = paste(col_vector[3]), shape="diamond")  
    graph <- visNetwork::visGroups(graph, groupname = "1Lpath", color = paste(col_vector[2]), shape="diamond")  
    graph <- visNetwork::visOptions(graph, highlightNearest = TRUE, nodesIdSelection = TRUE)
    graph <- visNetwork::visLegend(graph, useGroups = T)
    visNetwork::visSave(graph, file = paste(out.file,".html", sep=""))
  }else{
    query <- as.data.frame(gsub("map", species, FL), stringsAsFactors = F)
    colnames(query)="path_ID"
    #create store objects
    interactions <- NULL
    nextlevel <- NULL
    exclude <- NULL
    for (j in 2:length(level)){
      exclude <- rbind(exclude, query)
      colnames(exclude) = "path_ID"
      exclude <- unique(exclude)
      for (k in 1:length(query[,1])){
        out_nodes <- suppressMessages(data.frame(do.call(rbind, (XML::getNodeSet((XML::xmlParse(KEGGREST::keggGet(query[k,1], "kgml"))), "//entry[@type='map']", fun=XML::xmlToList)))))
        lks <- as.data.frame(t(as.data.frame(out_nodes$.attrs)), stringsAsFactors = F)
        lks <- unique(as.data.frame(paste(lks$name), stringsAsFactors = F))
        colnames(lks) <- "path_ID"
        #lks <- as.data.frame(subset(lks, !(lks$path_ID %in% query$path_ID)), stringsAsFactors = F)
        nextlevel <- rbind(nextlevel, lks)
        lks$to <- lks$path_ID
        if (length(lks$to)>0) {
          lks$from <- query[k,1]
        }else{
          lks$from <- lks$path_ID
        }
        lks <- lks[grep(species, lks$to), ]
        lks <- lks[lks$from != lks$to,]
        interactions <- rbind(interactions, lks)
      } 
      nextlevel <- unique(nextlevel)
      #grab only path with species specified
      nextlevel <- nextlevel[grep(species, nextlevel$path_ID, ignore.case = TRUE), ]
      #grab path different from previous level
      nextlevel <- nextlevel[(!(nextlevel %in% exclude$path_ID))]
      #path excluded in the next levels analysis
      exclude <- rbind(exclude, data.frame(path_ID = nextlevel))
      #create new query
      query <- data.frame(path_ID=nextlevel)
    }
  }
  interactions$path_ID <-NULL
  interactions <- unique(interactions)
  interactions <- interactions[interactions$from != interactions$to,]
  dm <- NULL
  for (i in 1:length(1:levels)){
    dm <- rbind(dm, get(paste("the",i,"Lgenes",sep="")))
  }
  dm <- as.data.frame(cbind(dm$gene_name, dm$path_ID), stringsAsFactors = F)
  colnames(dm) <- c("from", "to")
  int <- as.data.frame(cbind(interactions$from, interactions$to), row.names = F)
  colnames(int) <- c("from", "to")
  dm <- as.data.frame(rbind(dm, int))
  #extract the interactions uncovered to generate the diagram label
  label <- data.frame(name = unique(c(paste(dm$from), paste(dm$to))))
  label$num <- c(1:length(label$name))
  #create the attributes of diagram objects
  label$type <- 1
  for (k in length(1:levels):1){
    add <- get(paste("the",k,"Lgenes",sep=""))
    label$type <- ifelse(label$name %in% add$gene_name, paste(k, "Lgenes", sep=""), label$type)
    label$type <- ifelse(label$name %in% add$path_ID, paste(k, "Lpath", sep=""), label$type)
  }
  label$shape <- 1
  for (k in 1:length(1:levels)){
    add <- get(paste("the",k,"Lgenes",sep=""))
    label$shape <- ifelse(label$name %in% add$gene_name, "circle", label$shape)
  }
  label$shape <- ifelse(label$shape == "1", "diamond", label$shape)
  label$color <- "black"
  label$type <- ifelse(label$type == "1", "path without results", label$type)
  label$type <- ifelse(label$shape == "circle", "gene", label$type)
  label <- label[order(label$type, decreasing = F),]
  label <- label[order(label$shape, decreasing = F),]
  label <- label[order(label$num, decreasing = F),]
  #create the diagram legend
  legenda <- data.frame(cbind(num = label$num, name = paste(label$name)), stringsAsFactors = F)
  legenda <- merge(legenda, KEGGpath, by.x = "name", by.y = "rn", all=T)
  legenda <- legenda[order(as.numeric(legenda$num)),]
  legenda <- subset(legenda, (!is.na(legenda$num)))
  legenda$name <- ifelse(is.na(legenda$V1), legenda$name, legenda$V1)
  #create diagram nodes
  ndf <- data.frame(
    id = legenda$num,
    label = paste(legenda$name),
    shape = label$shape,
    group = label$type)
  #extract the interactions uncovered 
  for (i in 1:length(label$num)){
    dm$from[paste(dm$from) == paste(label$name[i])] <- paste(label$num[i])
    dm$to[paste(dm$to) == paste(label$name[i])] <- paste(label$num[i])
  }
  #create the diagram edges
  edf <- data.frame(
    from = dm$from,
    to = dm$to)
  #legend color
  qual_col_pals = RColorBrewer::brewer.pal.info[RColorBrewer::brewer.pal.info$category == 'qual',]
  col_vector = unlist(mapply(RColorBrewer::brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  #create the diagram and export a html file
  graph <- visNetwork::visNetwork(ndf, edf, main = paste("PANEV visualization", out.file, sep=" "))  
  graph <- visNetwork::visGroups(graph, groupname = "gene", color = paste(col_vector[1]), shape="dot")  
  graph <- visNetwork::visGroups(graph, groupname = "path without results", color = paste(col_vector[3]), shape="diamond")  
  graph <- visNetwork::visGroups(graph, groupname = "1Lpath", color = paste(col_vector[2]), shape="diamond")  
  graph <- visNetwork::visGroups(graph, groupname = "2Lpath", color = paste(col_vector[4]), shape="diamond")  
  graph <- visNetwork::visGroups(graph, groupname = "3Lpath", color = paste(col_vector[5]), shape="diamond")  
  graph <- visNetwork::visGroups(graph, groupname = "4Lpath", color = paste(col_vector[6]), shape="diamond")  
  graph <- visNetwork::visGroups(graph, groupname = "5Lpath", color = paste(col_vector[7]), shape="diamond")  
  graph <- visNetwork::visGroups(graph, groupname = "6Lpath", color = paste(col_vector[8]), shape="diamond")  
  graph <- visNetwork::visGroups(graph, groupname = "7Lpath", color = paste(col_vector[9]), shape="diamond")  
  graph <- visNetwork::visGroups(graph, groupname = "8Lpath", color = paste(col_vector[10]), shape="diamond")  
  graph <- visNetwork::visGroups(graph, groupname = "9Lpath", color = paste(col_vector[11]), shape="diamond")  
  graph <- visNetwork::visGroups(graph, groupname = "10Lpath", color = paste(col_vector[12]), shape="diamond")  
  graph <- visNetwork::visGroups(graph, groupname = "11Lpath", color = paste(col_vector[13]), shape="diamond")  
  graph <- visNetwork::visGroups(graph, groupname = "12Lpath", color = paste(col_vector[14]), shape="diamond")  
  graph <- visNetwork::visGroups(graph, groupname = "13Lpath", color = paste(col_vector[15]), shape="diamond")  
  graph <- visNetwork::visGroups(graph, groupname = "14Lpath", color = paste(col_vector[16]), shape="diamond")  
  graph <- visNetwork::visGroups(graph, groupname = "15Lpath", color = paste(col_vector[17]), shape="diamond")  
  graph <- visNetwork::visGroups(graph, groupname = "16Lpath", color = paste(col_vector[18]), shape="diamond")  
  graph <- visNetwork::visGroups(graph, groupname = "17Lpath", color = paste(col_vector[19]), shape="diamond")  
  graph <- visNetwork::visGroups(graph, groupname = "18Lpath", color = paste(col_vector[20]), shape="diamond")  
  graph <- visNetwork::visGroups(graph, groupname = "19Lpath", color = paste(col_vector[21]), shape="diamond")  
  graph <- visNetwork::visGroups(graph, groupname = "20Lpath", color = paste(col_vector[22]), shape="diamond")  
  graph <- visNetwork::visOptions(graph, highlightNearest = TRUE, nodesIdSelection = TRUE)
  graph <- visNetwork::visLegend(graph, useGroups = T)
  visNetwork::visSave(graph, file = paste(out.file,".html", sep=""))
  #clear the workspace
  panev.cleanEnvir(KEGGpath)
  panev.cleanEnvir(KEGGgenes)
  panev.cleanEnvir(store)
  setwd(wd)
  cat("\n")
  cat("Well done! Diagram visualization was created and exported. \n")
}