#Script: panev.uploadPANEVscriptoutput.R
#License: GPLv3 or later
#Modification date: 2019-05-13
#Written by: Valentino Palombo
#Contact: valentino.palombo@gmail.com
#Description: Upload the output(s) of PANEV script before perform the diagram visualization

panev.uploadPANEVscriptoutput <- function(
  levels = 1
)
{
  #upload the output file of PANEV script
  for (i in 1:length(1:levels)){
    store <<- paste("the",i, "Lgenes", sep="")
    assign(store, read.table(paste(i, "Lgenes.txt", sep=""), header = T, stringsAsFactors = F), envir = .GlobalEnv)
  }
}