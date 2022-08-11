

dargs <- list(i = "", o = "", b = FALSE)
Param <- R.utils::commandArgs(defaults = dargs, trailingOnly = FALSE, always = list(file =""), asValues = TRUE)



#Help message
if( any(names(Param) %in% c("h","help") ) ){
  cat("PAMLST a blast based tool for determining sequence type of pseudomonas aeruginosa genomes\n
Available commands:\n")
  cat("-i Input directory path Defaults to ./input/ in PAMLST folder if not specified
-o Output directory path (defaults to . if not specified)
-b (If you already previously performed the blast alignment for debugging purposes, defaults to FALSE)
-p (To enable installation of missing packages required to run script)
-u update MLST database and alleles
-h or -help Display this message\n")
  q()
  
}


#Set correct path for some functions
DIR2 <- dirname((Param$file))


if( !dir.exists( (DIR2) ) ){
  DIR3 <- ("./")
} else {
  DIR3 <- paste0( (DIR2),"/")
}

if( !dir.exists(DIR3) ){
  stop(paste0("something went wrong with path ",DIR3))
}



OsaDBdir <- paste0(DIR3,"Functions/OSAdb.fasta")
source(paste0(DIR3,"Functions/StupidMLST.R"))
dbDir <- paste0(DIR3,"db")








packages <- c("readr","utils")


#check for packages 

for(i in 1:length(packages)){
  a <- suppressWarnings(suppressMessages( require( packages[i], character.only= TRUE) ))
  
  if( !a & (! any(names(Param) %in% "p") )){
    stop(cat("Package", packages[i], "is missing, install package or run PAMLST.R with flag -p to enable installation of packages\n"))
  } else if( !a & any(names(Param) %in% "p") ){
    install.packages(packages[i])
    suppressWarnings(suppressMessages( require( packages[i]) ))
  }
  
}


if( any(names(Param) %in% c("u") ) ){
  print("updating MLST alleles and profile")
  source(paste0(DIR3,"Functions/DownloadMLST.R"))

  DownloadMLST(dbDir)

  #Concatenate alleles to a single multifasta file
  if(.Platform$OS.type == "unix"){
    system(paste0("cd ", dbDir, "; cat acsA.fas aroE.fas guaA.fas mutL.fas nuoD.fas ppsA.fas trpE.fas > MLSTdb.fasta;" ))

  } else if(.Platform$OS.type == "windows"){
    shell(paste0("cd ",dbDir, " & copy acsA.fas+aroE.fas+guaA.fas+mutL.fas+nuoD.fas+ppsA.fas+trpE.fas MLSTdb.fasta" ))
  } else {
    print("uh")
  }

  source(paste0(DIR3,"Functions/MLSTdbbuilder.R"))
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
  dInputquery <- paste0(DIR3,"Input/")
  

  
}
#Trailing slash to input
if( ! grepl("/$",dInputquery) ){
  dInputquery <- paste0(dInputquery, "/")
  
}
if( identical(list.files(dInputquery ,pattern = paste0(c("fasta$","fna$","fsa$","fas$","fa$","seq$"),collapse="|")), character(0)) ){
  stop("Input folder has no valid files, input files must have file extension: .fasta, .fna, .fsa, .fas, .fa, or .seq")
  # print(dOutput)
}


if( !identical( Param$o, "" ) ){
  
  dOutput <- Param$o
  if( !dir.exists(dOutput) ){
    dir.create(dOutput)
  }
  
} else {
  dOutput <- paste0(dirname(dInputquery),"/PAMLSToutput/" )
  
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

