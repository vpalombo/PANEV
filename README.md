Pathway Network Visualizer (PaNeV)
==================================

Package: PANEV

Type: R Package

Version: 1.0

Date: 2019-05-13

Author: Valentino Palombo, Marco Milanesi, Gabriella Sferra, Stefano Capomaccio, Sandy Sgorlon, Mariasilvia D'Andrea

Maintainer: Valentino Palombo <valentino.palombo@gmail.com>

Description: PANEV is an R package set for pathway-based network gene visualization. Based on information available on KEGG, it maps and visualizes genes within a network of upstream and downstream-connected pathways (from 1 to *n* levels). The graph helps to interpret functional profiles of cluster of genes.

License: Artistic-2.0

Abstract
========

Nowadays, the increasing in data generated from genomic and transcriptomic analyses requires strategies to solve the challenge of data mining. In this sense, the use of tools for pathway analysis and/or functional enrichment is *de facto* standard for the post analytic interrogation. In particular, when complex phenomenon is considered, create a network of multiple interconnected pathways of interest could be useful to investigate the underlying biology and, ultimately, identify functional candidate genes affecting the trait under investigation.

This manual describes the Pathways Network Visualizer (PANEV) R package. The idea is to create a 'functional' network between a set of main pathways of interest (first level - FL or 1L) and *n* levels of upstream/downstream interconnected pathways, to identify which genes, from a given list of genes of interest, are involved. The network refers to connected pathways and genes partecipating to sequential steps and thus presumed to affect the phenotype/condition of interest.

More in details, PANEV is based on information available on KEGG databases. The network is created based on the number (*n*) of levels required for the investigation. Starting from the FL (1L) chosen by the user, PANEV uncovers the connections between 1L and upstream/downstream pathways (from 1 to *n* levels). Then, once the backbone of pathways is created, PANEV highlights genes, among a list of provided ones, belonging to the multiple investigated levels and provides a pathways/genes-based network visualization. The highlighted genes may be considered good candidates for the trait/condition of interest, because supported by functional evidence. Overall, PANEV combined with classical functional enrichment analysis and/or pathway analysis can facilitate the interpretation of complex biological processes.

Since the package relies on KEGG, PANEV returns meaningful results only for genes with functional information available. The package is able to analyse genomic or transcriptomic outcomes for all species annotated in KEGG. The access to KEGG repositories has specific copyright conditions (<https://www.kegg.jp/kegg/legal.html>). PANEV uses KEGGREST package (<https://bioconductor.org/packages/release/bioc/html/KEGGREST.html>) functions to download individual pathway graphs and data files through API or HTTP access, which is freely available for academic and non-commercial uses.

Installation
============

This manual focuses on PANEV v.1.0. PANEV requires the pre-installation of ‘devtools’ package from CRAN (<https://cran.r-project.org/package=devtools>). Once satisfied the requirement, the R scripts and package described in this manual can be downloaded and installed from github link (<https://github.com/vpalombo/panev>) using:

    #install devtools
    install.packages("devtools") #if needed

    #install PANEV
    library("devtools")
    install_github("vpalombo/PANEV")

To load the package into the R environment (R ≥ 3.5.0) type:

    library("PANEV")

PANEV functions and workflow
============================

PANEV functions could be divided in two different steps: data preparation and data analyses.

Since PANEV interrogates biomaRt and/or KEGG databases, a **working internet connection** is required to properly run PANEV functions.

<br> <img src="https://github.com/vpalombo/PANEV/tree/master/vignettes/images/fig1.jpg" alt="Fig.1: The general architecture of the workflow of PANEV package and schematic illustration of main functions. The yellow rectangles represent the PANEV functions. The green circles represent the input data lists, in particular gene or pathway lists. The red diamonds represent the output from data preparation PANEV functions. The blue rectangles represent the final PANEV outcomes." width="450"> <br>

Example files
=============

Trial data files used in the help pages can be created with the command:

``` r
# Copy the example data in the current working directory
panev.example()
```

The files comprise:

-   for genomic dataset:
    -   *entrez\_genelist.txt* and *ensembl\_genelist.txt*. Two example files for genomic data preparation with entrez ID and ensembl ID annotation, respectively;
    -   *data.txt*. An example of gene list properly formatted ready to run PANEV main functions.
-   For transcriptomic dataset:
    -   *entrez\_expr\_genelist.txt* and *ensembl\_expr\_genelist.txt*. Two example files for transcriptomic data preparation with relative FC and p-value in entrez ID and ensembl ID annotation, respectively;
    -   *exprdata.txt*. An example of differentally expressed genes (DEG) list properly formatted and ready to run PANEV main functions;
    -   *expr\_listPath.txt*. An example file containing the list of main pathway of interest (FL or 1L) with relative expression estimated scores \[i.e. flux values obtained with Dinamic Impact Approach (Bionaz et al., 2012)\].

Input files used for the validation in reference publication can be created using the command:

``` r
# Copy the example data in the current working directory
panev.example(type="validation")
```

In this case, the files comprise:

-   *genelist\_annotated\_qiu2014.txt*. The gene list from a publicly available dataset containing 166 out of 171 'novel' genes discovered by a gene-base association study on human type 1 diabetes mellitus - T1DM (Qiu et al., 2014), properly formatted and ready to run PANEV functions.

-   *genelist\_enrichment\_qui2014.txt*. The gene list containing all the 452 genes discovered by Qiu et al. (2014).

-   *genelist\_expr\_Levy2012.txt*. Part of the previous gene list with relative FC and p-values, from a publicly available expression results obtained by Levy et al. (2012) and used by Qiu et al. (2014);

-   *pathlist\_expr\_Levy2012.txt*. The list of FL pathway.

Data preparation
================

The preparation of a properly formatted **list of genes** and **main pathways of interest** (FL or 1L) is the first step before running main PANEV functions.

Gene list preparation
---------------------

Since PANEV relies on KEGG databases information, the **entrez gene identifiers** are needed.

To endhance user experience, two specific functions, *panev.dataPreparation* and *panev.exprdataPreparation*, are provided to retrieve the correct ID annotation for genomic and transcriptomic data, respectively. Those functions are based on biomaRt R package that uses the *ensembl.org* website resource, that allows up-to-date identifiers switching.

The ensembl correct organism code is needed to run properly data preparation functions: these codes are obtainable with *panev.biomartSpecies*.

To create a list of all available organisms for biomaRt annotation:

``` r
# Create a list of all available organisms for biomaRt annotation
list <- panev.biomartSpecies(string = NULL)
#   The list of all available species for gene id convertion and data preparation was created! 
#   Remember to use the correct organism code for relative PANEV functions. 
```

To obtain the list of available code for a specific organism:

``` r
# Look for a specific organism matching a search string for biomaRt annotation
list <- panev.biomartSpecies(string = "cow")
#   The list of available species matched your string was created! 
#   Remember to use the correct organism code for relative PANEV functions.
list
#            organism_code        description version
#   1 btaurus_gene_ensembl Cow genes (ARS-UCD1.2)  ARS-UCD1.2
biomart.species.bos <- as.character(list[1,1])

list <- panev.biomartSpecies(string = "pig")
#   The list of available species matched your string was created! 
#   Remember to use the correct organism code for relative PANEV functions.
list
#                organism_code                           description     version
#   1     caperea_gene_ensembl Brazilian guinea pig genes (CavAp1.0)    CavAp1.0
#   2  cporcellus_gene_ensembl          Guinea Pig genes (Cavpor3.0)   Cavpor3.0
#   3 mnemestrina_gene_ensembl   Pig-tailed macaque genes (Mnem_1.0)    Mnem_1.0
#   4     sscrofa_gene_ensembl               Pig genes (Sscrofa11.1) Sscrofa11.1
biomart.species.sus <- as.character(list[4,1]) 
```

The use of data preparation functions is recommended but not mandatory. In fact, it must be emphasized that their correct performance depends on the availability to biomaRt data access for the desired species and on its correct entrez and/or ensembl ID annotation. For this reason, we strongly suggest to double-check all-possible genes annotation.

### For genomic dataset analysis

For genomic data, PANEV requires a specific format dataset. The list of the genes of interest must be provided as '*.txt*' file, stored in working directory, containing three columns labelled as *ensembl\_gene\_id*, *entrezgene* and *external\_gene\_name*, respectively.

User can create a properly formatted dataset from a simple gene list using *panev.dataPreparation* function. This function will convert the single gene list from **ensembl** to **entrez** gene ID (or vice versa), and add at the same time the **gene symbol**.

It is recommended to run this function before using *panev.script*.

Preparation of formatted gene list from ensembl gene ID:

``` r
# Example of gene list from genomic dataset with ensembl id 
genelist <- read.table("ensembl_genelist.txt")
head(genelist)
#                     V1
#   1 ENSBTAG00000000039
#   2 ENSBTAG00000000040
#   3 ENSBTAG00000000042
#   4 ENSBTAG00000000044
#   5 ENSBTAG00000001521
#   6 ENSBTAG00000001522

# Prepare the dataset for panev.script function converting gene name from ensembl to entrez id
genelist.converted <- panev.dataPreparation(in.file = "ensembl_genelist.txt", 
                                          gene_id = "ensembl", 
                                          biomart.species = biomart.species.bos)
#   Input file imported! 
#   
#   BiomaRt species correct! 
#   
#   Gene id correct!  
#   
#   Convertion from ensembl ID to entrez ID ... 
#   DONE
#   n. 37 of 47 genes have corresponding gene in KEGG database.
#   Gene list exported!
head(genelist.converted)
#        ensembl_gene_id entrezgene external_gene_name
#   1 ENSBTAG00000000039     505662              SIRT7
#   2 ENSBTAG00000000040     515219               MAFG
#   3 ENSBTAG00000000042     539606              PYCR1
#   4 ENSBTAG00000000044     617922            MYADML2
#   5 ENSBTAG00000001521     616871              UQCRB
#   6 ENSBTAG00000001522     526138             MTERF3
```

Preparation of formatted gene list from entrez gene ID:

``` r
# Example of gene list from genomic dataset with entrez id 
genelist <- read.table("entrez_genelist.txt")
head(genelist)
#            V1
#   1 104564307
#   2 100313156
#   3 100297382
#   4 100125314
#   5    782391
#   6    781828

# Prepare the dataset for panev.script function converting gene name from entrez to ensembl id
genelist.converted <- panev.dataPreparation(in.file = "entrez_genelist.txt", 
                                          gene_id = "entrez", 
                                          biomart.species = biomart.species.bos)
#   Input file imported! 
#   
#   BiomaRt species correct! 
#   
#   Gene id correct! 
#   
#   Convertion from entrez ID to ensembl ID ... 
#   DONE
#   n. 34 of 34 genes have corresponding gene in KEGG database.
#   Gene list exported!
head(genelist.converted)
#     entrezgene    ensembl_gene_id external_gene_name
#   1  100125314 ENSBTAG00000046173              ALG12
#   2  100297382 ENSBTAG00000021566               PAX2
#   3  100313156 ENSBTAG00000044643       bta-mir-2346
#   4  104564307 ENSBTAG00000039738             TMIGD3
#   5     280880 ENSBTAG00000013898                NPB
#   6     281073 ENSBTAG00000007591               CHUK
```

### For transcriptomic dataset analysis

For transcriptomic data, PANEV requires a specific format dataset. The gene list, containing the **differentially expressed genes** (DEG) of interest, must be provided as '*.txt*' file, stored in working directory, with five columns labelled as *ensembl\_gene\_id*, *entrezgene*, *external\_gene\_name*, *FC* (fold change) and *pvalue*.

User can create a PANEV formatted dataset using *panev.exprdataPreparation* function. This function will convert the DEG list with relative **FC** and **p-value** from **ensembl** to **entrez** gene ID (or vice versa), and add at the same time **gene symbols**.

It is recommended to run this function before using *panev.exprscript*.

Preparation of formatted DEG list from ensembl gene ID:

``` r
# Example of gene list from transcriptomic dataset with ensembl id 
genelist <- read.table("ensembl_expr_genelist.txt", header=TRUE)
head(genelist)
#        ensembl_gene_id           FC      pvalue
#   1 ENSSSCG00000000042     3.377911 0.000712043
#   2 ENSSSCG00000011956  -1.00335002 0.093817937
#   3 ENSSSCG00000027351  -1.90736265 0.012180974
#   4 ENSSSCG00000026626  1.360528698 0.083137308
#   5 ENSSSCG00000011621  2.824720966 0.000442323
#   6 ENSSSCG00000000042    3.3779108 0.000712043

# Prepare the dataset for panev.script function converting gene name from ensembl to entrez id
genelist.expr.converted <- panev.exprdataPreparation(in.file = "ensembl_expr_genelist.txt", 
                                                   gene_id = "ensembl", 
                                                   biomart.species = biomart.species.sus)
#   Input file imported! 
#   
#   BiomaRt species correct! 
#   
#   Gene id correct! 
#   
#   Convertion from ensembl ID to entrez ID ... 
#   DONE
#   n. 3708 of your 3774 genes have corresponding gene in KEGG database. 
#   Gene list exported!
head(genelist.expr.converted)
#        ensembl_gene_id          FC         pvalue   entrezgene external_gene_name
#   1 ENSSSCG00000000002   -4.992506     0.00002740   100520739               GTSE1
#   2 ENSSSCG00000000006 3.737968329     0.01211484      397239               PPARA
#   3 ENSSSCG00000000037 1.877643135    0.000200513   100523962             ARFGAP3
#   4 ENSSSCG00000000038  -1.201390429  0.052669253   100524254              CYB5R3
#   5 ENSSSCG00000000046  2.121964954   0.005759356      397687             CYP2D25
#   6 ENSSSCG00000000058  -1.426235589  0.001286279   100155420               SNU13
```

Preparation of formatted DEG list from entrez gene ID:

``` r
# Example of gene list from transcriptomic dataset with entrez id 
genelist <- read.table("entrez_expr_genelist.txt", header=TRUE)
head(genelist)
#     entrezgene         FC        pvalue
#   1     397239   3.737968   0.012114844
#   2  100523962   1.877643   0.000200513
#   3  100524254  -1.201390   0.052669253
#   4     397687   2.121965   0.005759356
#   5  100155420  -1.426236   0.001286279
#   6  100154200  -1.492203   0.000286301

# Prepare the dataset for panev.script function converting gene name from entrez to ensembl id
genelist.expr.converted <- panev.exprdataPreparation(in.file = "entrez_expr_genelist.txt", 
                                                   gene_id = "entrez", 
                                                   biomart.species = biomart.species.sus)
#   Input file imported! 
#   
#   BiomaRt species correct! 
#   
#   Gene id correct! 
#   
#   Convertion from ensembl ID to entrez ID ... 
#   DONE
#   n. 3531 of your 3531 genes have corresponding gene in KEGG database. 
#   Gene list exported!
head(genelist.expr.converted)
#     entrezgene         FC      pvalue      ensembl_gene_id  external_gene_name
#   1     387287  64.592010  0.000002350  ENSSSCG00000003601              HCRTR1
#   2     396564   1.248483  0.036568915  ENSSSCG00000012817              FUNDC2
#   3     396568   4.089142  0.006518320  ENSSSCG00000011802                KNG1
#   4     396574  -9.125993  0.006317556  ENSSSCG00000013387             
#   5     396575   1.679256  0.002075938  ENSSSCG00000000272                 SP1
#   6     396590   3.067011  0.003039118  ENSSSCG00000008857               MSMO1 
```

Pathways list preparation
-------------------------

The list of main pathways of interest (first level pathways - FL or 1L), coded as **KEGG path id**, is mandatory to run the main PANEV functions. The list of available pathways on KEGG is obtainable by *panev.pathList* function, that returns a dataframe containing two columns labelled as *path\_description* and *path\_ID*, respectively. This function helps to retrive the KEGG path ID if unknown.

To obtain the list of all available pathways in KEGG:

``` r
# Create a list of all available pathways for PANEV
list <- panev.pathList(string = NULL)
#   The list of all available pathways was created! 
#   Remember to use the correct path Id(s) for relative PANEV functions. 
```

To obtain the list of available pathway(s) matching a search string:

``` r
# Look for a specific pathway(s) for PANEV, matching your search string
list <- panev.pathList(string = "lipid")
#   The list of pathway(s), matched your string, was created! 
#   Remember to use the correct path Id(s) for relative PANEV functions.
list
#                                               path_description       path_ID
#   1                                    Glycerolipid metabolism path:map00561
#   2                             Glycerophospholipid metabolism path:map00564
#   3                                     Ether lipid metabolism path:map00565
#   4                                    Sphingolipid metabolism path:map00600
#   5 Glycosphingolipid biosynthesis - lacto and neolacto series path:map00601
#   6 Glycosphingolipid biosynthesis - globo and isoglobo series path:map00603
#   7            Glycosphingolipid biosynthesis - ganglio series path:map00604
#   8                             Sphingolipid signaling pathway path:map04071
#   9                                    Antidyslipidemic agents path:map07052
```

### For genomic dataset analysis

In case of genomic dataset analysis, the FL (1L) must be provided as a **vector of PANEV pathways identifiers** (*path\_ID*) obtained with *panev.pathList*.

``` r
# Create a vector of main pathways of interest (FL)
FL.gene <- c("path:map00061", "path:map00062", "path:map00071", "path:map00072")
```

### For transcriptomic dataset analysis

In case of transcriptomic dataset analysis, the (1L) FL must be provided as a **'*.txt*' file**, stored in working directory, with the **list of *path\_ID***, obtained with *panev.pathList*, and **relative expression estimated scores**, obtained by gene set enrichment analysis \[e.g. flux value (Bionaz et al., 2012)\].

``` r
# Example of main pathways of interest (FL) for transcriptomic data analysis
FL.expr <- read.table("expr_listPath.txt", header=TRUE)
FL.expr
#           path_ID     value
#   1 path:map00010 324.22879
#   2 path:map00020 -21.31287
#   3 path:map00071 385.73774
```

Species code detection
----------------------

The **KEGG** correct **organism code** is needed to run properly data analyses functions. The code is obtainable with *panev.speciesCode*.

To create a list of all available organisms in KEGG:

``` r
# Create a list of all available organisms in KEGG
list <- panev.speciesCode(string = NULL)
#   The list of all species available for PANEV analysis was created! 
#   Remember to use the correct organims code for relative PANEV functions.

# Number of available species in KEGG to date:
length(list$panev_code)
# [1] 5913

head(list)
# panev_code                                            species
#      hsa                                Homo sapiens (human)
#      ptr                        Pan troglodytes (chimpanzee)
#      pps                               Pan paniscus (bonobo)
#      ggo   Gorilla gorilla gorilla (western lowland gorilla)
#      pon                   Pongo abelii (Sumatran orangutan)
#      nle Nomascus leucogenys (northern white-cheeked gibbon)
```

To obtain the list of available code for a specific organism:

``` r
# Look for a specific organism in KEGG matching a search string
list <- panev.speciesCode(string = "bos")
#   The list of available species matched your string was created! 
#   Remember to use the correct organism code for relative PANEV functions.

head(list)
#                       species panev_code
#   1          Bos taurus (cow)      bta
#   2      Bos mutus (wild yak)      bom
#   3 Bos indicus (zebu cattle)      biu
#   4        Malassezia globosa      mgl
#   5      Bosea sp. PAMC 26642      bop
#   6           Bosea sp. RAC05      bos
KEGG.species.bos <- as.character(list[1,2])

list <- panev.speciesCode(string = "pig")
#   The list of available species matched your string was created! 
#   Remember to use the correct organism code for relative PANEV functions.

head(list)
#                         species panev_code
#   1            Sus scrofa (pig)      ssc
#   2 Columba livia (rock pigeon)      clv
#   3  Cajanus cajan (pigeon pea)     ccaj
#   4        Pigmentiphaga sp. H8      pig
#   5         Desulfovibrio piger      dpg
#   6         Salipiger profundus     tpro
KEGG.species.sus <- as.character(list[1,2]) 
```

Data analyses
=============

PANEV analysis and visualization
--------------------------------

PANEV automates the process of functional classification of genes of interest, detecting the genes inside a **functional network** of pathways. The network, created from selected main pathways (FL) and upstream/downstream dependent ones, helps to interpret functional profiles of genes, underlying complex biological processes.

### Genomic dataset

For genomic data, functional network is created, taking into account **upstream and downstream pathways** connected with FL pathways. The network is formed by pathways involved in sequential steps, connected at more levels (from **1** to **n**), as required by the user. After that, PANEV highlights the genes inside the network and provides a pathways/genes-based network visualization. The genes highlighted may be considered good candidates for the trait/condition of interest.

The *panev.script* function requires:

-   properly formatted gene list,

-   FL (1L), as KEGG pathways identifiers,

-   species code,

-   number of levels.

Please see the above "Data preparation" section to understand how to prepare the correct input files.

``` r
### Run if it is necessary ###
# Copy the example data file 'data.txt' in your current working directory
# panev.example()

# Head of input file "data.txt"
#        ensembl_gene_id entrezgene external_gene_name
#   1 ENSBTAG00000000039     505662              SIRT7
#   2 ENSBTAG00000000040     515219               MAFG
#   3 ENSBTAG00000000042     539606              PYCR1
#   4 ENSBTAG00000000044     617922            MYADML2
#   5 ENSBTAG00000001521     616871              UQCRB
#   6 ENSBTAG00000001522     526138             MTERF3

# Perform PANEV 
panev.script(in.file = "data.txt", 
           out.file = "FA", 
           species = KEGG.species.bos, 
           FL = FL.gene, 
           levels = 2)
#   Input file imported! 
#   
#   Gene list specified... and correct! 
#   
#   Species code specified... and correct! 
#   
#   Pathway(s) is specified... and correct! 
#   
#   Prerequisite check passed!
#   
#   PANEV is running ... 
#    Please wait... It could be a while depending on the number of pathways and levels required! 
#   
#   PANEV analysis completed and relative '.txt' files exported! 
#   
#   Preparing PANEV diagram visualization! 
#    Please wait... It could be a while depending on the number of pathways and levels required! 
#   
#   Well done! Diagram visualization was created and exported.       
```

A new folder is created and files resulting from the analysis are saved inside.

The function generates *n* text files, with *n* equal to the number of levels required for the investigation. Each single *.txt* file is a table containing the genes and the related pathways for a specific level

``` r
#Summary of PANEV results at FL level obtained with the example dataset
genes.1L <- read.table("PANEV_RESULTS_FA/1Lgenes.txt", header = TRUE)
genes.1L
#            ensemblgene entrezgene gene_name        path_description       path_ID
#   1 ENSBTAG00000015980     281152      FASN Fatty acid biosynthesis path:bta00061
#   2 ENSBTAG00000015178     505355      ECI2  Fatty acid degradation path:bta00071
```

Along with the **tabular format**, the function gets also the **diagram visualization** of PANEV results, saved in an interactive '*.html*' file. The diagram allows zooming on all content for an optimal readability.

<br> <img src="https://github.com/vpalombo/PANEV/tree/master/vignettes/images/fig2.jpg" alt="Fig.2: An example of the 'html' file with the network-based visualization of PANEV results. The green circles represent the  functional candidate genes inside the network of pathways generated. The violet diamonds represent the first level (FL) pathway(s), directly connected to the trait of interest and showing candidate gene(s). The yellow diamonds represent the second level of pathways connected with FL pathways and showing candidate gene(s). The orange diamonds represent the pathways investigated without any candidate gene. The diagram allows identifying the relationships among genes and pathways for each investigated level." width="450"> <br>

The '*html*' diagram is **interactive and interrogable**. You can select nodes by ID/label (e.g *FASN* gene).

<br> <img src="https://github.com/vpalombo/PANEV/tree/master/vignettes/images/fig3.jpg" alt="Fig.3: An example of node selection of PANEV network-based visualization result." width="450"> <br>

### Transcriptomic dataset

For transcriptomic data, PANEV takes into account any possible connection **only among the FL pathways** and a list of differentially expressed genes (DEG).

The *panev.script* function requires:

-   properly formatted DEG list,

-   properly formatted FL (1L) list,

-   species code,

-   p-value cut-off (for filtering subset of genes in the DEG list).

Please see the above "Data preparation" section to understand how to prepare the correct input files.

``` r
### Run if it is necessary ###
# Copy the example data files "exprdata.txt" and 'expr_listPath.txt' in your current working directory
# panev.example()

# Head of "exprdata.txt" file
#        ensembl_gene_id         FC        pvalue  entrezgene   external_gene_name
#   1 ENSSSCG00000000002  -4.992506   0.00002740    100520739                GTSE1
#   2 ENSSSCG00000000006   3.737968  0.012114844       397239                PPARA
#   3 ENSSSCG00000000037   1.877643  0.000200513    100523962              ARFGAP3
#   4 ENSSSCG00000000038  -1.201390  0.052669253    100524254               CYB5R3
#   5 ENSSSCG00000000046   2.121965  0.005759356       397687              CYP2D25
#   6 ENSSSCG00000000058  -1.426236  0.001286279    100155420                SNU13

# Head of "expr_listPath.txt" file
#           path_ID     value
#   1 path:map00010 324.22879
#   2 path:map00020 -21.31287
#   3 path:map00071 385.73774

# Perform PANEV on transcriptomic dataset
panev.exprscript(in.file = "exprdata.txt", 
               path.file = "expr_listPath.txt", 
               out.file = "expression_data", 
               species = KEGG.species.sus, 
               pvalue = 0.05)
#   Input file imported! 
#   
#   Pathway input file imported! 
#   
#   Gene list specified... and correct! 
#   
#   Your path list colnames are correct! 
#   
#   Species code specified... and correct! 
#   
#   Pathway(s) is specified... and correct! 
#   
#   Prerequisite check passed!
#   
#   PANEV is running ... 
#    Please wait... It could be a while depending on the number of pathways required! 
#   n. 3132 of 3774 genes passed the p-value filtering.
#   
#   Well done! Diagram visualization was created and exported. 
```

A new folder is created and files resulting from the analysis are saved inside.

The function generates a diagram visualization, taking into account any possible connections among the FL pathways and DEGs. The diagram helps to interpret the results obtained from gene expression experiments showing the nodes (i.e. genes and pathways) coloured according to their FCs and pathway expression estimated scores \[e.g. flux values (Bionaz et al., 2012)\], respectively.

The classification into upregulated/downregulated genes/pathways is done as function of top FC/estimated score values.

| gene/pathway classification |                with FC/estimated score value                |
|:---------------------------:|:-----------------------------------------------------------:|
|     low up/downregulated    |         &lt;25% of top up/downregulated gene/pathway        |
|  moderate up/downregulated  | &gt;= 25% and &lt; 50% of top up/downregulated gene/pathway |
|    high up/downregulated    | &gt;= 50% and &lt; 75% of top up/downregulated gene/pathway |
|   strong up/downregulated   |        &gt;= 75% of top up/downregulated gene/pathway       |

The '*html*' diagram is interactive and allows zooming on all content for an optimal readability.

<br> <img src="https://github.com/vpalombo/PANEV/tree/master/vignettes/images/fig4.jpg" alt="Fig.4: An example of the 'html' file with the network-based visualization of PANEV result considering an expression dataset. The circles represent the genes coloured based on their fold change (FC) values. The diamonds represent the pathways of interest coloured based on their expression estimated scores [i.e. flux values obtained with Dinamic Impact Approach (Bionaz et al., 2012)]. The diagram shows the relationships among genes and pathways and allows identifying functionally related entities with possibly coordinated expression changes." width="450"> <br>

The '*html*' diagram is also interrogable. You can select nodes by ID/label (e.g. *ALDH2* gene).

<br> <img src="https://github.com/vpalombo/PANEV/tree/master/vignettes/images/fig5.jpg" alt="Fig.5: Example of node selection of PANEV network-based visualization result obtained on transcriptomic dataset." width="450"> <br>

Enrichment analysis
-------------------

PANEV can perform also an enrichment analysis for each KEGG term (i.e. pathways) based on **hypergeometric test** (one-sided Fisher exact test) as described by Simoes and Emmert-Streib (2012).

The results are a series of '*.txt*' files with specific enrichment analysis results and with general descriptive information about gene and pathway occurrences.

More in details:

-   \*&lt;out.file&gt;\_enrichment.txt\* contains the enrichment analysis result. A p-value is calculated for each pathway to estimate the probability of significant over-representation.

-   \*&lt;out.file&gt;\_GxP.txt\* contains the list of pathways with gene occurrences.

-   \*&lt;out.file&gt;\_PxG.txt\* contains the list of genes with pathway occurrences.

Those files are useful to explore the list of interesting pathways which might be possible candidates for further investigation, particularly when no previous restrictive biological assumptions on trait/condition of interest are available.

``` r
### Run if it is necessary ###
# Copy the example data file 'data.txt' in your current working directory
# panev.example()

# Head of input file "data.txt"
#        ensembl_gene_id entrezgene external_gene_name
#   1 ENSBTAG00000000039     505662              SIRT7
#   2 ENSBTAG00000000040     515219               MAFG
#   3 ENSBTAG00000000042     539606              PYCR1
#   4 ENSBTAG00000000044     617922            MYADML2
#   5 ENSBTAG00000001521     616871              UQCRB
#   6 ENSBTAG00000001522     526138             MTERF3

# Perform the enrichment analysis
panev.stats.enrichment(in.file = "data.txt", 
                     out.file = "enrichment_FA", 
                     species = KEGG.species.bos)
#   Input file is imported! 
#   
#   Gene list specified... and correct! 
#   
#   Species code specified... and correct! 
#   
#   Enrichment analysis started ... 
#   and results exported! 
#   
#   Gene per pathway(s) table created and exported! 
#   
#   Pathway per gene(s) table created and exported!

# <out.file>_enrichment.txt file
FA_enrich <- read.table("enrichment_FA_enrichment.txt", header = T)
head(FA_enrich)
#        pathway_ID n_genes all_genes       pvalue       padj                           pathway_name
#   1 path:bta01100       8      1466  0.002440842  0.6758655                     Metabolic pathways
#   2 path:bta04920       2        72  0.007663162  0.6758655        Adipocytokine signaling pathway
#   3 path:bta05212       2        74  0.008080410  0.6758655                      Pancreatic cancer
#   4 path:bta04662       2        75  0.008292828  0.6758655      B cell receptor signaling pathway
#   5 path:bta00440       1         6  0.011094948  0.7233906 Phosphonate and phosphinate metabolism
#   6 path:bta04152       2       123  0.021290522  1.0000000                 AMPK signaling pathway

# <out.file>_GxP.txt file
FA_GxP <- read.table("enrichment_FA_GxP.txt", header = T)
head(FA_GxP)
#     n_genes               pathway_name    pathway_ID
#   1       8         Metabolic pathways path:bta01100
#   2       2     MAPK signaling pathway path:bta04010
#   3       2      Ras signaling pathway path:bta04014
#   4       2     FoxO signaling pathway path:bta04068
#   5       2 PI3K-Akt signaling pathway path:bta04151
#   6       2     AMPK signaling pathway path:bta04152

# <out.file>_PxG.txt file
FA_PxG <- read.table("enrichment_FA_PxG.txt", header = T)
head(FA_PxG)
#     n_pathways entrez_gene_id    ensembl_gene_id gene_symbol
#   1         42         281073 ENSBTAG00000007591        CHUK
#   2         23         619066 ENSBTAG00000022927        RAC3
#   3         12         369023 ENSBTAG00000016253       G6PC3
#   4          8         616871 ENSBTAG00000001521       UQCRB
#   5          5         281152 ENSBTAG00000015980        FASN
#   6          3         510274 ENSBTAG00000001868       PCYT2
```

Validation in publication
=========================

A publicly available data set from a study on human type 1 diabetes mellitus - T1DM (Qiu et al., 2014) was used as testing set for PANEV in our publication.

In the reference study, the authors carried out a gene-based genome-wide association analysis and identified 452 significant genes. Among these genes, 171 were newly identified for type 1 diabetes mellitus, since ignored in previously studies or in literature. The authors further supported the significance of 53 out 171 discovered genes with replication studies and differential expression studies.

In particular, the authors reported four non-HLA genes (*RASIP1*, *STRN4*, *BCAR1* and *MYL2*) and three *HLA* genes (*FYN*, *HLA-J* and *PPP1R11*) as validated by both replication and differential expression studies and represent the main result discussed by the authors.

We performed PANEV considering the list of 171 newly genes, to verify the possible contribution of PANEV visualization for candidate genes identification and more broadly for high-throughput data interpretation.

Genomic dataset analysis
------------------------

'*The Type I diabetes mellitus*', '*Insulin resistance*' and '*AGE-RAGE signaling pathway in diabetic complications*' were chosen as FL pathways, since clearly connected to the trait of interest (Greenbaum, 2002; Ramasamy et al., 2005).

Considering the complexity of the trait, the analysis was performed up to the third level (Field et al., 1997).

After data preparation, only 5 out of 171 genes had no entrez gene identifier. This list was used for running the PANEV function.

``` r
### Run if it is necessary ###
# Copy the example files used as validation set in the publication in the current working directory
# panev.example(type="validation")

# Parameters
in.file="genelist_annotated_qiu2014.txt"
# Head of in.file
#   external_gene_name entrezgene ensembl_gene_id
# 1              ADAD1     132612 ENSG00000164113
# 2              ASCL2        430 ENSG00000183734
# 3             ATF7IP      55729 ENSG00000171681
# 4               BAK1        578 ENSG00000030110
# 5              BCAR1       9564 ENSG00000050820
# 6             BCL2A1        597 ENSG00000140379

#Look for the specie code matching the search string 
list <- panev.speciesCode(string = "homo")
species=as.character(list[1,2]) # hsa 

FL = c("path:map04940", "path:map04931", "path:map04933")
levels = 3


# Run the PANEV function
panev.script(in.file = in.file, 
           out.file = "validation", 
           species = species, 
           FL = FL, 
           levels = levels)
#   Input file imported! 
#   
#   Gene list specified... and correct! 
#   
#   Species code specified... and correct! 
#   
#   Pathway(s) is specified... and correct! 
#   
#   Prerequisite check passed!
#   
#   PANEV is running ... 
#    Please wait... It could be a while depending on the number of pathways and levels required! 
#   
#   PANEV analysis completed and relative '.txt' files exported! 
#   
#   Preparing PANEV diagram visualization! 
#    Please wait... It could be a while depending on the number of pathways and levels required! 
#   
#   Well done! Diagram visualization was created and exported. 
```

Overall PANEV results obtained from validation dataset are in line with reference study outcomes (Qiu et al., 2014), confirming the effectiveness of our visualization approach.

In particular, 4 out of 7 genes validated both replication and differential expression studies (Qui et al. 2014) were highlighted by PANEV: *PTPN11*, *BCAR1*, *MYL2* and *FYN*. The other three genes (*RASIP1*, *STRN4* and *HLA-J*) were not detected by PANEV, since, although present in KEGG databases, they were not yet assigned to any pathway.

Along with these genes, PANEV highlighted also other interesting genes (*ITPR3*, *BAK1*, *IL10*, *HMGB1* and *MICA*) not discussed by Qui et al. (2014), since validated only by the differential expression study or only by the replication study.

``` r
#Summary of PANEV results considering first level (FL or 1L) pathways
genes.1L <- read.table("PANEV_RESULTS_validation/1Lgenes.txt", header = TRUE)
genes.1L
#     ensemblgene       entrezgene  gene_name        path_description  path_ID
# 1   ENSG00000179295       5781    PTPN11          Insulin resistance  path:hsa04931
```

``` r
#Summary of PANEV results considering second level (2L) pathways
genes.2L <- read.table("PANEV_RESULTS_validation/2Lgenes.txt", header = TRUE)
genes.2L
#     ensemblgene       entrezgene  gene_name        path_description  path_ID
# 1 ENSG00000123374       1017      CDK2        PI3K-Akt signaling pathway path:hsa04151
# 2 ENSG00000096433       3710     ITPR3                         Apoptosis path:hsa04210
# 3 ENSG00000030110        578      BAK1                         Apoptosis path:hsa04210
# 4 ENSG00000140379        597    BCL2A1                         Apoptosis path:hsa04210
# 5 ENSG00000010810       2534       FYN T cell receptor signaling pathway path:hsa04660
# 6 ENSG00000136634       3586      IL10 T cell receptor signaling pathway path:hsa04660
# 7 ENSG00000179295       5781    PTPN11   Adipocytokine signaling pathway path:hsa04920
# 8 ENSG00000101665       4092    SMAD7   TGF-beta signaling pathway path:hsa04350
```

``` r
#Summary of PANEV results considering third level (3L) pathways
genes.1L <- read.table("PANEV_RESULTS_validation/1Lgenes.txt", header = TRUE)
genes.1L
#      ensemblgene  entrezgene  gene_name                            path_description       path_ID
# 1  ENSG00000096433       3710     ITPR3                   Calcium signaling pathway path:hsa04020
# 2  ENSG00000050820       9564     BCAR1                 Chemokine signaling pathway path:hsa04062
# 3  ENSG00000140379        597    BCL2A1                NF-kappa B signaling pathway path:hsa04064
# 4  ENSG00000136634       3586      IL10                      FoxO signaling pathway path:hsa04068
# 5  ENSG00000123374       1017      CDK2                      FoxO signaling pathway path:hsa04068
# 6  ENSG00000096433       3710     ITPR3       Phosphatidylinositol signaling system path:hsa04070
# 7  ENSG00000123374       1017      CDK2                                  Cell cycle path:hsa04110
# 8  ENSG00000123374       1017      CDK2                       p53 signaling pathway path:hsa04115
# 9  ENSG00000189403       3146     HMGB1                          Autophagy - animal path:hsa04140
# 10 ENSG00000030110        578      BAK1 Protein processing in endoplasmic reticulum path:hsa04141
# 11 ENSG00000010810       2534       FYN                              Focal adhesion path:hsa04510
# 12 ENSG00000050820       9564     BCAR1                              Focal adhesion path:hsa04510
# 13 ENSG00000111245       4633      MYL2                              Focal adhesion path:hsa04510
# 14 ENSG00000099866       8174   MADCAM1              Cell adhesion molecules (CAMs) path:hsa04514
# 15 ENSG00000179295       5781    PTPN11                  Jak-STAT signaling pathway path:hsa04630
# 16 ENSG00000136634       3586      IL10                  Jak-STAT signaling pathway path:hsa04630
# 17 ENSG00000138378       6775     STAT4                  Jak-STAT signaling pathway path:hsa04630
# 18 ENSG00000204520  100507436      MICA   Natural killer cell mediated cytotoxicity path:hsa04650
# 19 ENSG00000010810       2534       FYN   Natural killer cell mediated cytotoxicity path:hsa04650
# 20 ENSG00000179295       5781    PTPN11   Natural killer cell mediated cytotoxicity path:hsa04650
# 21 ENSG00000111245       4633      MYL2            Regulation of actin cytoskeleton path:hsa04810
# 22 ENSG00000050820       9564     BCAR1            Regulation of actin cytoskeleton path:hsa04810
```

PANEV highlighted also other genes not discussed in the reference study (Qiu et al., 2014) but reported in literature as being associated to the susceptibility to T1DM disease. In particular, *CDK2* (Kim et al., 2017), *MADCAM1* (Phillips et al., 2005), *STAT4* (Bi et al., 2013), *SMAD7* (Chen et al., 2011) and *BCL2A1* (Beyan et al. 2010).

It is worth to note that, except for *CDK2*, all genes mentioned above refer to researches conducted before the reference study (Qiu at al., 2014).

Simultaneously, it is must be observed that some genes were not highlighted by PANEV, because:

-   fell inside no-investigated pathway (*BRAP*, *FUT2*, *GNS*, *HIPK1*, *NUPR1*, *OR2B3*, *HIST1H4E*, *HIST1H2BF*, *OR2B3*, *OR2B6*, *OR2J2*, *OR5V1*, and *SULT1A1* genes). About 8% of our investigated genes.

-   although present in KEGG databases, were not yet assigned to any pathway. About 44% of our investigated genes.

These drawbacks clearly represent the main PANEV limitation.

Enrichment analysis
-------------------

Accordly to reference study (Qui et al., 2014), we performed the functional annotation clustering analysis of 452 identified T1DM genes, using the PANEV enrichment analysis function.

``` r
### Run if it is necessary ###
# Copy the example files used as validation set in the publication in the current working directory
# panev.example(type="validation")

# Parameters
in.file="genelist_enrichment_qui2014.txt"
# Head of in.file
#   external_gene_name entrezgene ensembl_gene_id
# 1             OLFML3      56944   not_available
# 2              HIPK1     204851   not_available
# 3               IL10       3586   not_available
# 4               NSL1      25936   not_available
# 5             FAM46B     115572   not_available
# 6               LHX9      56956   not_available

#Look for the specie code matching the search string 
list <- panev.speciesCode(string = "homo")
species=as.character(list[1,2]) # hsa 

# Run the PANEV enrichment function
panev.stats.enrichment(in.file = in.file, 
                     out.file = "validation", 
                     species=species)
#   Input file is imported! 
#   
#   Gene list specified... and correct! 
#   
#   Species code specified... and correct! 
#   
#   Enrichment analysis started ... 
#   and results exported! 
#   
#   Gene per pathway(s) table created and exported! 
#   
#   Pathway per gene(s) table created and exported!
```

The results showed as the genes tend to be over-represented in immune diseases and immune system pathways, according to Qiu et al. (2014).

``` r
#Summary of PANEV enrichment result
enrichment.result <- read.table("validation_enrichment.txt", header = TRUE)
head(enrichment.result)
#      pathway_ID n_genes all_genes       pvalue         padj                        pathway_name
# 1 path:hsa05322      26       133 1.338552e-22 4.403835e-20        Systemic lupus erythematosus
# 2 path:hsa05330      13        38 3.566961e-15 5.867652e-13                 Allograft rejection
# 3 path:hsa04612      16        77 1.306151e-14 1.432412e-12 Antigen processing and presentation
# 4 path:hsa05320      14        53 1.763583e-14 1.450547e-12          Autoimmune thyroid disease
# 5 path:hsa04940      13        43 2.264972e-14 1.490352e-12            Type I diabetes mellitus
# 6 path:hsa05332      12        41 3.705738e-13 2.031980e-11           Graft-versus-host disease
```

Transciptomic dataset analysis
------------------------------

By way of example, to explore the possibility to visualize gene expression values, we used in the publication the FC values obtained by Levy et al. (2012) study and considered by Qiu et al. (2014) as reference study for differential expression validation.

Since the authors did not provide pathway estimated scores, we substituted those measures with gene occurrence values for each pathway of interest, obtained by PANEV enrichment analysis.

``` r
### Run if it is necessary ###
# Copy the example data files "genelist_expr_Levy2012.txt" and 'pathlist_expr_Levy2012.txt' in your current working directory
# panev.example(type="validation")

# Parameters
in.file="genelist_expr_Levy2012.txt"
#
# Head of "genelist_expr_Levy2012" file
#  ensembl_gene_id     external_gene_name   entrezgene   FC     pvalue 
#  ENSG00000204252     HLA-DOA              3111        -0.535  0.000
#  ENSG00000239457     HLA-DOB              3112        -0.017  0.955
#  ENSG00000168384     HLA-DPA1             3113        -0.491  0.025
#  ENSG00000206239     HLA-DQA1             3117        -0.846  0.000
#  ENSG00000206237     HLA-DQB1             3119        -0.467  0.007
#  ENSG00000204592     HLA-E                3133         0.298  0.014
#
path.file = "pathlist_expr_Levy2012.txt"
# Head of "pathlist_expr_Levy2012" file
# path_ID         value
# path:map04514    16
# path:map04940    14
# path:map04151     6
# path:map04210     5
# path:map04630     5
# path:map04010     4

#Look for the specie code matching the search string 
list <- panev.speciesCode(string = "homo")
species=as.character(list[1,2]) # hsa

# Perform PANEV on transcriptomic dataset
panev.exprscript(in.file = in.file, 
               path.file = path.file, 
               out.file = "expression_data_validation", 
               species = species, 
               pvalue = 0.05)
#   Input file imported! 
#   
#   Pathway input file imported! 
#   
#   Gene list specified... and correct! 
#   
#   Your path list colnames are correct! 
#   
#   Species code specified... and correct! 
#   
#   Pathway(s) is specified... and correct! 
#   
#   Prerequisite check passed!
#   
#   PANEV is running ... 
#    Please wait... It could be a while depending on the number of pathways required! 
#   n. 18 of 30 genes passed the p-value filtering.
#   
#   Well done! Diagram visualization was created and exported. 
```

This PANEV visualization helped to interpret the results obtained from gene expression experiments by showing the nodes (i.e. genes) coloured according to gene FC values. The diagram showed the relationships among genes and pathways and allowed us to identify functionally related entities with possibly coordinated expression changes.

Conclusion
==========

PANEV represents an interesting visualization approach for reducing the complexity of the challenge of high-throughput data mining and for candidate genes identification. PANEV providing an integrated summary among pathways and genes of interest, facilitating the data interpretation.

The contribution of PANEV tool could be significant not only for well-annotated species (i.e. Homo sapiens, Mus musculus) but also for all the organisms available in KEGG databases. Although KEGG is a popular and constantly updated database for biological network information, the lack or incomplete information on KEGG could represent the main PANEV disadvantage. The effectiveness of PANEV analysis in terms of result coherency was confirmed by the validation study. In particular, PANEV produces timesaving advantages, pointing users to genes that are biologically involved with the investigated trait.

References
==========

Beyan,H. et al. (2010) Monocyte gene-expression profiles associated with childhood-onset type 1 diabetes and disease risk: a study of identical twins. Diabetes, 59, 1751–1755. <doi:10.2337/db09-1433>.

Bi,C. et al. (2013) Association study of STAT4 polymorphisms and type 1 diabetes in Northeastern Chinese Han population. Tissue Antigens, 81, 137–140. <doi:10.1111/tan.12057>.

Bionaz,M. et al. (2012). A novel dynamic impact approach (DIA) for functional analysis of time-course omics studies: validation using the bovine mammary transcriptome. PLos One 7:e32455. <doi:10.1371/journal.pone.0032455>.

Chen HY, Huang XR, Wang W, Li JH, Heuchel RL, Chung ACK, et al. The protective role of Smad7 in diabetic kidney disease: mechanism and therapeutic potential. Diabetes. 2011;60: 590-601. <doi:10.2337/db10-0403>

Field,L.L. at al. (1997). Unravelling a complex trait: the genetics of insulin-dependent diabetes mellitus. Clin. Investig. Med. Med. Clin. Exp. 20:41–49.

Greenbaum, C.J. 2002. Insulin resistance in type 1 diabetes. Diabetes Metab. Res. Rev. 18:192–200. <doi:10.1002/dmrr.291>.

Kim,S.Y. et al. (2017) Loss of cyclin dependent kinase 2 in the pancreas links primary β-cell dysfunction to progressive depletion of β-cell mass and diabetes. J. Biol. Chem., <doi:10.1074/jbc.M116.754077>.

Levy, H., X. Wang, M. Kaldunski, S. Jia, J. Kramer, S.J. Pavletich, M. Reske, T. Gessel, M. Yassai, M.W. Quasney, M.K. Dahmer, J. Gorski, and M.J. Hessner. 2012. Transcriptional signatures as a disease-specific and predictive inflammatory biomarker for type 1 diabetes. Genes Immun. 13:593-604. <doi:10.1038/gene.2012>.

Phillips,J.M. et al. (2005) MAdCAM-1 is needed for diabetes development mediated by the T cell clone, BDC-2·5. Immunology, 116, 525–531. <doi:10.1111/j.1365-2567.2005.02254.x>.

Qiu,Y.-H. et al. (2014) Identification of novel risk genes associated with type 1 diabetes mellitus using a genome-wide gene-based association analysis. J. Diabetes Investig., 5, 649–656. <doi:10.1111/jdi.12228>.

Ramasamy, R. et al. (2005). Advanced glycation end products and RAGE: a common thread in aging, diabetes, neurodegeneration, and inflammation. Glycobiology 15:16R-28R. <doi:10.1093/glycob/cwi053>.

Shi,A. et al. (2016) Genetic variants in vitamin D signaling pathways and risk of gestational diabetes mellitus. Oncotarget, 7, 67788–67795. <doi:10.18632/oncotarget.11984>.

Simoes, R. de M., and F. Emmert-Streib. 2012. Bagging statistical network inference from large-scale gene expression data. PLos One 7:e33624. <doi:10.1371/journal.pone.0033624>.
