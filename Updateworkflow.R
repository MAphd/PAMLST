Usedir <- dirname(rstudioapi::getSourceEditorContext()$path)

#Updating alleles and profile:

#Download alleles
download.file("https://pubmlst.org/bigsdb?db=pubmlst_paeruginosa_seqdef&page=downloadAlleles&locus=acsA", paste0(Usedir,"/db/acsA.fas"))
download.file("https://pubmlst.org/bigsdb?db=pubmlst_paeruginosa_seqdef&page=downloadAlleles&locus=aroE", paste0(Usedir,"/db/aroE.fas"))
download.file("https://pubmlst.org/bigsdb?db=pubmlst_paeruginosa_seqdef&page=downloadAlleles&locus=guaA", paste0(Usedir,"/db/guaA.fas"))
download.file("https://pubmlst.org/bigsdb?db=pubmlst_paeruginosa_seqdef&page=downloadAlleles&locus=mutL", paste0(Usedir,"/db/mutL.fas"))
download.file("https://pubmlst.org/bigsdb?db=pubmlst_paeruginosa_seqdef&page=downloadAlleles&locus=nuoD", paste0(Usedir,"/db/nuoD.fas"))
download.file("https://pubmlst.org/bigsdb?db=pubmlst_paeruginosa_seqdef&page=downloadAlleles&locus=ppsA", paste0(Usedir,"/db/ppsA.fas"))
download.file("https://pubmlst.org/bigsdb?db=pubmlst_paeruginosa_seqdef&page=downloadAlleles&locus=trpE", paste0(Usedir,"/db/trpE.fas"))

#Download profile
download.file("https://pubmlst.org/bigsdb?db=pubmlst_paeruginosa_seqdef&page=downloadProfiles&scheme_id=1", paste0(Usedir,"/db/paeruginosamlstprofile.txt"))


#concenate MLST alleles


shell(paste0("cd ",Usedir, "/db/ & copy acsA.fas+aroE.fas+guaA.fas+mutL.fas+nuoD.fas+ppsA.fas+trpE.fas MLSTdb.fasta" ))

#Run MLSTdbbuilder 
source(paste0(Usedir,"/Functions/MLSTdbbuilder.R"))
