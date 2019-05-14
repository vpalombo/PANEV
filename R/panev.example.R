#' @title Copy the example files in the Working Directory
#' @description The function copies the example files necessary to test the package in the current working directory.
#' @usage panev.example(type=NULL)
#' @param type If set as "validation", the files used in the paper are saved in the working directory. In all the other cases, the example files used in the help section are saved in the working directory (default = NULL).
#' @details The example files used in the help section are:
#' \itemize{
#' \item{\emph{ensembl_genelist.txt}: this is the input file to test \code{\link[PANEV]{panev.dataPreparation}}. Please check function help to have information about the file format. }
#' \item{\emph{entrez_genelist.txt}:  this is the input file to test \code{\link[PANEV]{panev.dataPreparation}}. Please check function help to have information about the file format. }
#' \item{\emph{data.txt}: this is the input file to test \code{\link[PANEV]{panev.script}} and \code{\link[PANEV]{panev.stats.enrichment}}. Please check function help to have information about the file format. }
#' \item{\emph{ensembl_expr_genelist.txt}: this is the input file to test \code{\link[PANEV]{panev.exprdataPreparation}}. Please check function help to have information about the file format. }
#' \item{\emph{entrez_expr_genelist.txt}: it is the example input file to test \code{\link[PANEV]{panev.exprdataPreparation}}. Please check function help to have information about the file format. }
#' \item{\emph{expr_listPath.txt}: this is the input file to test \code{\link[PANEV]{panev.exprscript}}. Please check function help to have information about the file format. }
#' \item{\emph{exprdata.txt}: this is the input file to test \code{\link[PANEV]{panev.exprscript}} and \code{\link[PANEV]{panev.stats.enrichment}}. Please check functions help to have information about the file format. }
#' }
#' 
#' The files used in the publication to validate PANEV R package are: 
#' \itemize{
#' \item{\emph{genelist_annotated_qiu2014.txt}: this is the input file used in the publication to test \code{\link[PANEV]{panev.script}}. The gene list is from a publicly available data set from a study on human type 1 diabetes mellitus - T1DM (Qiu et al., 2014). }
#' \item{\emph{genelist_expr_Levy2012.txt}: this is the input file for "in.file" argument used in the publication to test \code{\link[PANEV]{panev.exprscript}}. The gene list with relative FC values is from  expression results obtained by Levy et al. (2012). }
#' \item{\emph{pathlist_expr_Levy2012.txt}: this is the input file for "path.file" argument used in the publication to test \code{\link[PANEV]{panev.exprscript}}. The pathway list with relative expression score values (in this case subsitute by number of occurrence) is from expression results obtained by Levy et al. (2012). }
#' }
#' In the vignette is reported how to replicate the analyses presented in the publication.  
#' @author Marco Milanesi (\email{marco.milanesi.mm@gmail.com})
#' @references Levy, H., X. Wang, M. Kaldunski, S. Jia, J. Kramer, S.J. Pavletich, M. Reske, T. Gessel, M. Yassai, M.W. Quasney, M.K. Dahmer, J. Gorski, and M.J. Hessner. 2012. Transcriptional signatures as a disease-specific and predictive inflammatory biomarker for type 1 diabetes. Genes Immun. 13:593-604. doi:10.1038/gene.2012
#' @references Qiu, Y.-H., F.-Y. Deng, M.-J. Li, and S.-F. Lei. 2014. Identification of novel risk genes associated with type 1 diabetes mellitus using a genome-wide gene-based association analysis. J. Diabetes Investig. 5:649-656. doi:10.1111/jdi.12228.
#' @examples ##### EXAMPLES CODE #####
#' #Copy the example files used in the help sections in the current working directory
#' panev.example()
#' 
#' #Copy the example files used as validation set in the publication in the current working directory
#' panev.example(type="validation")


#Script: panev.example.R
#License: GPLv3 or later
#Modification date: 2019-05-13
#Written by: Marco Milanesi
#Contact: marco.milanesi.mm@gmail.com
#Description: Copy the example files in the wd

panev.example<-function(type=NULL){
  #Copy files to current directory
  if (!(is.null(type))){
    if (type=="validation"){
      file.copy(system.file("data","genelist_annotated_qiu2014.txt", package = "PANEV"), "genelist_annotated_qiu2014.txt")
      file.copy(system.file("data","genelist_expr_Levy2012.txt", package = "PANEV"), "genelist_expr_Levy2012.txt")
      file.copy(system.file("data","pathlist_expr_Levy2012.txt", package = "PANEV"), "pathlist_expr_Levy2012.txt")
      file.copy(system.file("data","genelist_enrichment_qui2014.txt", package = "PANEV"), "genelist_enrichment_qui2014.txt")
    }else{
      file.copy(system.file("data","data.txt", package = "PANEV"), "data.txt")
      file.copy(system.file("data","ensembl_expr_genelist.txt", package = "PANEV"), "ensembl_expr_genelist.txt")
      file.copy(system.file("data","ensembl_genelist.txt", package = "PANEV"), "ensembl_genelist.txt")
      file.copy(system.file("data","entrez_expr_genelist.txt", package = "PANEV"), "entrez_expr_genelist.txt")
      file.copy(system.file("data","entrez_genelist.txt", package = "PANEV"), "entrez_genelist.txt")
      file.copy(system.file("data","expr_listPath.txt", package = "PANEV"), "expr_listPath.txt")
      file.copy(system.file("data","exprdata.txt", package = "PANEV"), "exprdata.txt") 
    } 
  }else{
    file.copy(system.file("data","data.txt", package = "PANEV"), "data.txt")
    file.copy(system.file("data","ensembl_expr_genelist.txt", package = "PANEV"), "ensembl_expr_genelist.txt")
    file.copy(system.file("data","ensembl_genelist.txt", package = "PANEV"), "ensembl_genelist.txt")
    file.copy(system.file("data","entrez_expr_genelist.txt", package = "PANEV"), "entrez_expr_genelist.txt")
    file.copy(system.file("data","entrez_genelist.txt", package = "PANEV"), "entrez_genelist.txt")
    file.copy(system.file("data","expr_listPath.txt", package = "PANEV"), "expr_listPath.txt")
    file.copy(system.file("data","exprdata.txt", package = "PANEV"), "exprdata.txt") 
  }
}