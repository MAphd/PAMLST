
#### MLST 

StupidMLST <- function(Input, Output, dbDir, Blasted) { 

#Determine the sequence type of a given genome using BLASTn
library(readr)
  

MLSTdb <- paste0(dbDir,"/MLSTdb.fasta")
MLSTdbtype <- paste0(dbDir,"/paeruginosamlstprofile.txt")
# Output <- dOutput
Alleledir <- dbDir
# Input <- dInput

  
  
DST <- data.frame(read_delim(MLSTdbtype, 
                                            "\t", escape_double = FALSE, trim_ws = TRUE, col_types = cols()))



# Check output folder exists, otherwise create it (dir.exists, dir.create)
# Ensure all paths end with a forwardslash , in general probably a good idea to test if all dirs are valid 
# 

# library(readr)


#Use folder names as ST for now

#So build a dataframe which lists the genome names, the associated serotype (PAst), and the sequence type

# grep("contigs", list.dirs(Input,recursive = FALSE))


# Inputs <- list.dirs(Input, recursive = FALSE)[c(grep("contigs", list.dirs(Input,recursive = FALSE)))]
# 
# #For now just 1
# Input <- Inputs[2]
# Input <- sub("//","/",Input)
# 
# Input <- paste0(Input,"/")

Seqnms <- list.files(Input)

#We could go trhough the pain of identifying filenames and valid filenames,, but we wont for now

# for(i in 1:length(Inputs))

if(length(Seqnms)>=1){ 
  print( paste("Succesfully loaded ",length(Seqnms)," sequences", sep=""))
} else {
  
  stop("No input files could be loaded")
}

if( !Blasted ){ 

Cmdline <- function(Seqfile){
  paste("blastn -query ", Input, Seqfile, " -subject ", MLSTdb, " -out ",
        Output, Seqfile,
        ".txt -outfmt ", "\"", "6 saccver pident nident mismatch gaps sstart send length evalue bitscore", "\"", " -max_target_seqs 5000",  sep="")
  
}
Cmds <- mapply(FUN = Cmdline, Seqfile = Seqnms) 

# 
# paste("blastn -query ", Input, Seqfile, " -subject ", MLSTdb, " -out ",
#       Output, Seqfile,
#       ".txt -outfmt ", "\"", "6 qaccver saccver pident sstart send length evalue bitscore", "\"",  sep="")


# system(Cmds[20])
#Run each blast colnames
sapply(Cmds, system)  

}


Data <- readRDS(paste0((Alleledir),"/MLSTdb.rds"))

DL <- list()
for( i in 1:length(ls(Data))){
  
  Temp <- eval(parse(text=(paste0("Data$",ls(Data)[i]))))
  
  DL[[ ls(Data)[i] ]] <- list()
  
  for( j in 1:length(Temp)){






    DL[[ ls(Data)[i] ]][[ names(Temp[[j]]) ]] <- length(unlist((strsplit(Temp[[j]],""))))


  }
  
  
}
rm(Data)


Blastnms <- paste(Seqnms,".txt",sep="")

Finaldata <- data.frame(matrix(nrow=length(Blastnms), ncol=9))
colnames(Finaldata) <- c("Seqname","ST","acsA","aroE","guaA","mutL","nuoD","ppsA","trpE")
# Finaldata <- Finaldata[-1,]
for(i in 1:length(Seqnms)){
 
    # i <- 1
  Data <- read.delim(paste(Output, Blastnms,sep="")[i], colClasses = "character", header=FALSE, col.names = c("Subjacc","Identity","Nident","Mismatch","Gaps","Sstart","Send","Length","Evalue","Bitscore"))
  #Compute coverage for each hit
  #For every Data$Subjacc

  # Refsize <- DL[[strsplit(Data$Subjacc[1],"_")[[1]][1]]][as.numeric(strsplit(Data$Subjacc[1],"_")[[1]][2])]

  # 100 - 100*(Refsize-Data)/Refsize

  Allrefsizes <- mapply( FUN = function(y) {

    DL[[strsplit(Data$Subjacc[y],"_")[[1]][1]]][(Data$Subjacc[y])][[1]]
    
    }, y=1:dim(Data)[1])

  #any(is.na(Allrefsizes))
  
  Data$Coverage <- (100 - 100*(as.numeric(Allrefsizes)-as.numeric(Data$Length) )/as.numeric(Allrefsizes)  )
  
  #Typing!
  
  #easy way first: 
  #Check if we have 7 hits, (and if these 7 are all unique)
  
  #Do we have 7 perfect hits? (0 mismatched, 0 gaps and 100% coverage (no more, no less))
  if( identical(dim((Data[Data$Mismatch==0&Data$Gaps==0&Data$Coverage==100,]))[1] , as.integer(7)) ){
    #are these 7 hits unique alleles (not just two perfect hits for 1 allele etc)
    a <- unlist(strsplit( sort(Data[Data$Mismatch==0&Data$Gaps==0&Data$Coverage==100,]$Subjacc) , "_" ))
    if(    identical(length(unique(a[seq(1,length(a),2)])),as.integer(7))  ){
      
      # a[seq(2,length(a),2)]
      
      seqST <- as.numeric(a[c(FALSE,TRUE)])
    } else {
      seqST <- as.numeric(c(1,0,0,0,0,0,0) ) # if 7 unique alleles cant be found
    }
    
    
  } else {
    # if we get more or less than 7 perfect hits 
    seqST <- as.numeric(c(0,1,0,0,0,0,0) )
    # Is there multiple hits for some alleles and more than 7 hits? last ditch effort to assign a ST
    if(  dim((Data[Data$Mismatch==0&Data$Gaps==0&Data$Coverage==100,]))[1] > as.integer(7) & 
         any(duplicated(unlist(strsplit(Data$Subjacc[Data$Mismatch==0&Data$Gaps==0&Data$Coverage==100],"_"))[c(TRUE,FALSE)]))
    ){
      
      #Index of duplicated allele
      b <- duplicated(unlist(strsplit(Data$Subjacc[Data$Mismatch==0&Data$Gaps==0&Data$Coverage==100],"_"))[c(TRUE,FALSE)])
      #if more than one duplicated allele, deem untypable (ie ST 9999)
      if( identical(sum(b),as.integer(1)) ){
        
        
        
        #check possible MLSTs 
        a <- strsplit(Data$Subjacc[Data$Mismatch==0&Data$Gaps==0&Data$Coverage==100],"_")
        #duplicated allele name
        a1 <- a[[which(b)]][[1]]
        
        #indexes of duplicated allele
        b1 <- which(unlist(strsplit(Data$Subjacc[Data$Mismatch==0&Data$Gaps==0&Data$Coverage==100],"_"))[c(TRUE,FALSE)]%in%a1)
        
        a2 <- unlist(strsplit( sort(Data[Data$Mismatch==0&Data$Gaps==0&Data$Coverage==100,]$Subjacc[-(b1[1])]) , "_" ))
        if(    identical(length(unique(a2[seq(1,length(a2),2)])),as.integer(7))  ){
          seqST1 <- as.numeric(a2[c(FALSE,TRUE)])
          
        }
        
        a3 <- unlist(strsplit( sort(Data[Data$Mismatch==0&Data$Gaps==0&Data$Coverage==100,]$Subjacc[-(b1[2])]) , "_" ))
        if(    identical(length(unique(a3[seq(1,length(a3),2)])),as.integer(7))  ){
          
          seqST2 <- as.numeric(a3[c(FALSE,TRUE)])
          
        }
        
        #Check if MLST exists for both allele combinations
        check1 <- DST$ST[DST$acsA==seqST1[1]&DST$aroE==seqST1[2]&DST$guaA==seqST1[3]&DST$mutL==seqST1[4]&DST$nuoD==seqST1[5]&DST$ppsA==seqST1[6]&DST$trpE==seqST1[7]]
        check2 <- DST$ST[DST$acsA==seqST2[1]&DST$aroE==seqST2[2]&DST$guaA==seqST2[3]&DST$mutL==seqST2[4]&DST$nuoD==seqST2[5]&DST$ppsA==seqST2[6]&DST$trpE==seqST2[7]]
        
        #Check if any MLSTs were assigned
        if( any((!c(identical(check1,numeric(0)),identical(check2,numeric(0))))) ){
          
          #If more than 1 MLSTs can be assigned, there are duplicate alleles and the genome is untypable so we don't proceed
          if( !identical(sum((!c(identical(check1,numeric(0)),identical(check2,numeric(0))))),as.integer(2)) ){
            
            if( identical(which(!c(identical(check1,numeric(0)),identical(check2,numeric(0)))),as.integer(1)) ){
              
              seqST <- as.numeric(a2[c(FALSE,TRUE)])
            } else if( identical(which(!c(identical(check1,numeric(0)),identical(check2,numeric(0)))),as.integer(2)) ){
              
              seqST <- as.numeric(a3[c(FALSE,TRUE)])
            }
            
            
          }
          
          
        }
        
        
        
        
        rm(a1);rm(a2);rm(a3);rm(b);rm(b1);rm(seqST1);rm(seqST2);rm(check1);rm(check2)
      }
    }
    
  }
  
  
  #lookup the typing scheme from paeruginosamlstprofile.txt and assign this to the sequence


  
  #make sure every allele has been assigned before trying this
  if( identical( length(seqST),as.integer(7))&(!any(seqST==0))  ){ 
  seqMLST <- DST$ST[DST$acsA==seqST[1]&DST$aroE==seqST[2]&DST$guaA==seqST[3]&DST$mutL==seqST[4]&DST$nuoD==seqST[5]&DST$ppsA==seqST[6]&DST$trpE==seqST[7]]
  # identical(as.numeric(DST[146,2:8]), as.numeric(seqST))
  if( identical(seqMLST, numeric(0)) ){
    #If we have a novel ST (not listed )
    seqMLST <- 9997
    
  }
  
  } else { 
    if( identical( seqST, c(0,1,0,0,0,0,0) ) ){
      seqMLST <- 9999
      
    } else{
      
      if( identical( seqST, c(1,0,0,0,0,0,0) ) ){
        seqMLST <- 9998
        
        
      }
      
      
    }
    
    
    
    }
  #Now we wanna start building our output
  
  Finaldata$Seqname[i] <- Blastnms[i]
  Finaldata$ST[i] <- seqMLST
  Finaldata$acsA[i] <- seqST[1]
  Finaldata$aroE[i] <- seqST[2]
  Finaldata$guaA[i] <- seqST[3]
  Finaldata$mutL[i] <- seqST[4]
  Finaldata$nuoD[i] <- seqST[5]
  Finaldata$ppsA[i] <- seqST[6]
  Finaldata$trpE[i] <- seqST[7]
  # Finaldata <- rbind(Finaldata, 
  # c(Blastnms[i], seqMLST, seqST
  # 
  # 
  #       
  # )
  #                                       
  # )
  
  
  
}

write.table(Finaldata, file=paste(Output,gsub("/","",strsplit(Input,dirname(Input))[[1]][2]),"_MLST.txt",sep=""),row.names=FALSE)

Finaldata

}
