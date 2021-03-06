% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/panev.example.R
\name{panev.example}
\alias{panev.example}
\title{Copy the example files in the Working Directory}
\usage{
panev.example(type=NULL)
}
\arguments{
\item{type}{If set as "validation", the files used in the paper are saved in the working directory. In all the other cases, the example files used in the help section are saved in the working directory (default = NULL).}
}
\description{
The function copies the example files necessary to test the package in the current working directory.
}
\details{
The example files used in the help section are:
\itemize{
\item{\emph{ensembl_genelist.txt}: this is the input file to test \code{\link[PANEV]{panev.dataPreparation}}. Please check function help to have information about the file format. }
\item{\emph{entrez_genelist.txt}:  this is the input file to test \code{\link[PANEV]{panev.dataPreparation}}. Please check function help to have information about the file format. }
\item{\emph{data.txt}: this is the input file to test \code{\link[PANEV]{panev.network}} and \code{\link[PANEV]{panev.stats.enrichment}}. Please check function help to have information about the file format. }
\item{\emph{ensembl_expr_genelist.txt}: this is the input file to test \code{\link[PANEV]{panev.exprdataPreparation}}. Please check function help to have information about the file format. }
\item{\emph{entrez_expr_genelist.txt}: it is the example input file to test \code{\link[PANEV]{panev.exprdataPreparation}}. Please check function help to have information about the file format. }
\item{\emph{expr_listPath.txt}: this is the input file to test \code{\link[PANEV]{panev.exprnetwork}}. Please check function help to have information about the file format. }
\item{\emph{exprdata.txt}: this is the input file to test \code{\link[PANEV]{panev.exprnetwork}} and \code{\link[PANEV]{panev.stats.enrichment}}. Please check functions help to have information about the file format. }
}

The files used in the publication to validate PANEV R package are: 
\itemize{
\item{\emph{genelist_annotated_qiu2014.txt}: this is the input file used in the publication to test \code{\link[PANEV]{panev.network}}. The gene list is from a publicly available data set from a study on human type 1 diabetes mellitus - T1DM (Qiu et al., 2014). }
\item{\emph{genelist_expr_Levy2012.txt}: this is the input file for "in.file" argument used in the publication to test \code{\link[PANEV]{panev.exprnetwork}}. The gene list with relative FC values is from  expression results obtained by Levy et al. (2012). }
\item{\emph{pathlist_expr_Levy2012.txt}: this is the input file for "path.file" argument used in the publication to test \code{\link[PANEV]{panev.exprnetwork}}. The pathway list with relative expression score values (in this case subsitute by number of occurrence) is from expression results obtained by Levy et al. (2012). }
}
In the vignette is reported how to replicate the analyses presented in the publication.
}
\examples{
##### EXAMPLES CODE #####
#Copy the example files used in the help sections in the current working directory
panev.example()

#Copy the example files used as validation set in the publication in the current working directory
panev.example(type="validation")
}
\references{
Levy, H., X. Wang, M. Kaldunski, S. Jia, J. Kramer, S.J. Pavletich, M. Reske, T. Gessel, M. Yassai, M.W. Quasney, M.K. Dahmer, J. Gorski, and M.J. Hessner. 2012. Transcriptional signatures as a disease-specific and predictive inflammatory biomarker for type 1 diabetes. Genes Immun. 13:593-604. doi:10.1038/gene.2012

Qiu, Y.-H., F.-Y. Deng, M.-J. Li, and S.-F. Lei. 2014. Identification of novel risk genes associated with type 1 diabetes mellitus using a genome-wide gene-based association analysis. J. Diabetes Investig. 5:649-656. doi:10.1111/jdi.12228.
}
\author{
Marco Milanesi (\email{marco.milanesi.mm@gmail.com})
}
