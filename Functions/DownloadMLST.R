DownloadMLST <- function(dbDir){


#Updating alleles and profile:

print("Downloading alleles")
download.file("https://pubmlst.org/bigsdb?db=pubmlst_paeruginosa_seqdef&page=downloadAlleles&locus=acsA", paste0(dbDir,"/acsA.fas"))
download.file("https://pubmlst.org/bigsdb?db=pubmlst_paeruginosa_seqdef&page=downloadAlleles&locus=aroE", paste0(dbDir,"/aroE.fas"))
download.file("https://pubmlst.org/bigsdb?db=pubmlst_paeruginosa_seqdef&page=downloadAlleles&locus=guaA", paste0(dbDir,"/guaA.fas"))
download.file("https://pubmlst.org/bigsdb?db=pubmlst_paeruginosa_seqdef&page=downloadAlleles&locus=mutL", paste0(dbDir,"/mutL.fas"))
download.file("https://pubmlst.org/bigsdb?db=pubmlst_paeruginosa_seqdef&page=downloadAlleles&locus=nuoD", paste0(dbDir,"/nuoD.fas"))
download.file("https://pubmlst.org/bigsdb?db=pubmlst_paeruginosa_seqdef&page=downloadAlleles&locus=ppsA", paste0(dbDir,"/ppsA.fas"))
download.file("https://pubmlst.org/bigsdb?db=pubmlst_paeruginosa_seqdef&page=downloadAlleles&locus=trpE", paste0(dbDir,"/trpE.fas"))

print("Downloading MLST profile")
download.file("https://pubmlst.org/bigsdb?db=pubmlst_paeruginosa_seqdef&page=downloadProfiles&scheme_id=1", paste0(dbDir,"/paeruginosamlstprofile.txt"))


#concenate MLST alleles



#Run MLSTdbbuilder 

}