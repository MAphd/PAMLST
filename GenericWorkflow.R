
library(readr)

#Load required functions

Usedir <- dirname(rstudioapi::getSourceEditorContext()$path)
#Attempt to define directory automatically



source(paste0(Usedir,"/Functions/StupidMLST.R"))

#Define variables for running PAst:
# output: output directory
# input: input directory
# dbDir: db directory 
# blasted: TRUE/FALSE - omits blasting step if set to TRUE

#Define these if not using original file structure
Inputdir <- paste0(Usedir,"/Input/") 
Outputdir <- paste0(Usedir,"/Output/")
dbDir <- paste0(Usedir,"/db/")
Blasted <- FALSE


StupidMLST(Inputdir, Outputdir, dbDir, Blasted)


#Output file will be saved to Outputdir with the suffix "_MLST.txt".



