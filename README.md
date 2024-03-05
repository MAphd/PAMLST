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

Sequence types 9997, 9998 and 9999 are internal error codes.
9997 is if the genome has a novel ST. 9998 is if 7 unique alleles can't be found. 9999 is if we get more or less than 7 perfect hits to the MLST alleles.

Required software:
* Rscript
* blastn 

Required R packages:
* readr

Has been tested using Rscript 3.6.3 for Ubuntu 20.04 with blastn 2.9.0

If you found this useful please cite https://doi.org/10.1099/mgen.0.000919
```
Anbo M, Jelsbak L. 2023. A bittersweet fate: detection of serotype switching in Pseudomonas aeruginosa. Microb Genomics 9:000919. 
```
