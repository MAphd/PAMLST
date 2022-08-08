
#Load functions, must be in same dir as this script
source("Functions/StupidMLST.R")







#Define these if not using original file structure
# Inputdir <- paste0(Usedir,"/Input/") 
# Outputdir <- paste0(Usedir,"/Output/")
dbDir <- "db"







dargs <- list(i = "", o = "", b = FALSE)
Param <- R.utils::commandArgs(defaults = dargs, trailingOnly = TRUE, asValues = TRUE)



#Help message
if( any(names(Param) %in% c("h","help") ) ){
  cat("PAMLST a blast based tool for determining sequence type of pseudomonas aeruginosa genomes\n
Available commands:\n")
  cat("-i Input directory path (Must be a file compatible with blastn with file extensions: .fasta, .fsa, .fas, .fna, .fa. or .seq). Defaults to Input/ if not specified
-o Output directory path (defaults to output/ if not specified)
-b (If you already previously performed the blast alignment for debugging purposes, defaults to FALSE)
-p (To enable installation of missing packages required to run script)
-u update MLST database and alleles
-h or -help Display this message\n")
  q()
  
}



packages <- c("readr","utils")


#check for packages 

for(i in 1:length(packages)){
  a <- suppressWarnings(suppressMessages( require( packages[i], character.only= TRUE) ))
  
  if( !a & (! any(names(Param) %in% "p") )){
    stop(cat("Package", packages[i], "is missing, install package or run PAMLST.R with flag -p to enable installation of packages\n"))
  } else if( !a & any(names(Param) %in% "p") ){
    install.packages(packages[i])
  }
  
}


if( any(names(Param) %in% c("u") ) ){
  print("updating MLST alleles and profile")
  source("Functions/DownloadMLST.R")

  DownloadMLST(dbDir)

  #Concatenate alleles to a single multifasta file
  if(.Platform$OS.type == "unix"){
    system(paste0("cd ", dbDir, "; cat acsA.fas aroE.fas guaA.fas mutL.fas nuoD.fas ppsA.fas trpE.fas > MLSTdb.fasta;" ))

  } else if(.Platform$OS.type == "windows"){
    shell(paste0("cd ",dbDir, " & copy acsA.fas+aroE.fas+guaA.fas+mutL.fas+nuoD.fas+ppsA.fas+trpE.fas MLSTdb.fasta" ))
  } else {
    print("uh")
  }

  source("Functions/MLSTdbbuilder.R")
  MLSTdbbuilder(dbDir)
  
  print("Successfully updated MLST alleles and profile")
  q()
  
}

system("blastn -version",intern=TRUE)





if( !identical( Param$i, "" ) ){
  
  dInputquery <- as.character(Param$i)
  if( !dir.exists(dInputquery) ){
    stop("Input invalid")
  }
  
} else {
  dInputquery <- "Input/"
  
  if( identical(list.files(dInputquery ,pattern = paste0(c("*fasta$","*fna$","*fsa$","*fas$","*fa$","*seq$"),collapse="|")), "as.character(0)") ){
    stop("Input folder has no valid files, input files must have file extension: .fasta, .fna, .fsa, .fas, .fa, or .seq")
    # print(dOutput)
  }
  
}
#Trailing slash to input
if( ! grepl("/$",dInputquery) ){
  dInputquery <- paste0(dInputquery, "/")
  
}




if( !identical( Param$o, "" ) ){
  
  dOutput <- Param$o
  if( !dir.exists(dOutput) ){
    dir.create(dOutput)
  }
  
} else {
  dOutput <- "Output/"
  
  if( !dir.exists(dOutput) ){
    dir.create(dOutput)
    # print(dOutput)
  }
  
}


if( !is.logical(Param$b) ){
  stop("Blasted must be either TRUE or FALSE or left blank (defaults to FALSE)")
  
} else {
  
  Blasted <- Param$b
}

StupidMLST(dInputquery, dOutput, dbDir, Blasted)

