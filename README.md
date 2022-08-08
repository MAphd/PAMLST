# PAMLST

A blast based tool for determining sequence type of pseudomonas aeruginosa genomes

Call from commandline using 

```
Rscript PAMLST.R 
```

Avaliable commands:

```
-i Input directory path (Must be a file compatible with blastn with file extensions: .fasta, .fsa, .fas, .fna, .fa. or .seq). Defaults to ./input/ if not specified
-o Output directory path (defaults to ./output/ if not specified)
-b (If you already previously performed the blast alignment for debugging purposes, defaults to FALSE) 
-p (To enable installation of missing packages required to run script)
-u to update MLST database and alleles
-h or -help Display this message
```

Required software:
* Rscript
* blastn 

Required R packages:
* readr

Has been tested using Rscript 3.6.3 for Ubuntu 20.04 with blastn 2.9.0

