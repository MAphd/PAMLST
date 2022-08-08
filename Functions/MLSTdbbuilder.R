MLSTdbbuilder <- function(Alleledir){


# library(readr)

# Usedir <- dirname(rstudioapi::getSourceEditorContext()$path)
# 
# Alleledir <- paste0(Usedir,"/db/")



DE <- new.env(hash = TRUE)
for( j in 1:length( list.files(Alleledir, "*.fas") )){
  
  # j <- 3
  Data <- suppressMessages(suppressWarnings(read_delim(paste0(Alleledir,"/", list.files(Alleledir, "*.fas")[j]),">",trim_ws = FALSE, col_names = FALSE, col_types = cols())))
  
  #Get indexes of where a sequence starts and stops
  d <- which(!is.na(Data$X2))
  
  
  #Create a list for each allele, listname is taken from the second coloumn, and the sequence is given to this.
  D <- list()
  #D2 <- new.env()
  for(i in 1:length(d)) {
    
    if ( i < length(d) ){
      b <- paste0(unlist(Data[(d[i]+1):(d[i+1]-1),1]),collapse="")
    } else {
      b <- paste0(unlist(Data[(d[i]+1):dim(Data)[1],1]),collapse="")
    }
    
    names(b) <- strsplit(as.character(Data[d[i],]),">")[[2]]
    
    D[[ strsplit(names(b),"_")[[1]][[2]] ]] <- b
    
    # D[[b]] <- strsplit(as.character(Data[d[i],]),">")[[2]]
    
    
  }
  
  # DL <- Kmerize(D,35,1)
  
  # names(DL)
  
  # paste0(strsplit(get(ls(DL[[1]])[1],DL[[1]]),"_")[[1]][1:2],sep="",collapse="_")
  
  
  
  
  # names(DL) <- mapply( FUN = function(nx) { paste0(strsplit(get(ls(DL[[nx]])[1],DL[[nx]]),"_")[[1]][1:2],sep="",collapse="_") }, nx=1:length(DL) )
  
  # DL <- as.environment(DL)
  
  #Append all alleles, ie append(DL[[1]],DL[[2]])
  #However, to do that, we need to remove any duplicates within a locus, (then later across all locus)
  # for(ij in 1:length(d)) {
  #   for(ji in 1:(length(d)-1) )
  #   which(duplicated())
  #   
  # 
  # }
  
  
  # for(ij in 1:3   ){
  #   
  #   ij <- 2
  #   if( ij+1 <= length(d)) {
  #   for(ji in (ij+1):length(d)) { 
  #     #Vector of duplicates between list ij and ji
  #     a <- which(duplicated( c(ls(DL[[ij]]),ls(DL[[ji]]) ) )) - length(ls(DL[[ij]]))
  #     #Then we want to remove the duplicates from list ji
  #     if( !identical(a,integer(0) ) ){
  #     DL[[ji]] <- DL[[ji]][-a]
  #     }
  #     
  #   }
  #   }
  #   
  # }
  
  # a <- which(duplicated(c(ls(DL[[1]]),ls(DL[[80]]))))- length(ls(DL[[1]]))
  # 
  # DL[[2]] <- DL[[2]][-a]
  # 
  # which( ls(DL[[2]]) == c(ls(DL[[1]]),ls(DL[[2]]))[357] )
  # 
  #  ( 357- length(ls(DL[[1]])) )
  
 # DE[[ strsplit(D[[1]],"_")[[1]][[1]] ]] <- D
  
  
  DE[[ strsplit(names(D[[1]]),"_")[[1]][1] ]] <- D
  
  
  
}

saveRDS(DE, file=paste0(Alleledir,"/MLSTdb.rds"))

}


# 
# print(paste0("Database made with ", length(DE), " different locus: ", (paste(ls(DE),sep=" ",collapse=" ")), collapse="" ))
# 
# print(paste("Saved database to", paste0(dirname(Alleledir), "/KmerMLST.rds" ))) 
# }