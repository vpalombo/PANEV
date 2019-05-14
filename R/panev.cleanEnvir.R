#Script: panev.cleanEnvir.R
#License: GPLv3 or later
#Modification date: 2019-05-13
#Written by: Valentino Palombo
#Contact: valentino.palombo@gmail.com
#Description: Function to clean the workspace

panev.cleanEnvir <- function(x) {rm(list=deparse(substitute(x)),envir=.GlobalEnv)}