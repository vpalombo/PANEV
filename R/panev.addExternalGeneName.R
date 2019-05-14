#Script: panev.addExternalGeneName.R
#License: GPLv3 or later
#Modification date: 2019-05-13
#Written by: Valentino Palombo
#Contact: valentino.palombo@gmail.com
#Description: Add gene id to empty gene symbol records.

panev.addExternalGeneName <- function(
  genelist=genelist
)
{
  #substitute the empty gene symbol record(s) with the ensembl or entrez gene id
  genelist[which(genelist$external_gene_name == ""), "external_gene_name"] <- genelist[which(genelist$external_gene_name == ""),"ensembl_gene_id"]
  genelist[which(genelist$external_gene_name == ""), "external_gene_name"] <- genelist[which(genelist$external_gene_name == ""),"entrezgene"]
  return(genelist)
}