## RECKONER - Read Error Corrector Based on KMC
### RECKONER
RECKONER is a tool for correction of Illumina reads. It is based
on the error correction algorithm BLESS and utilizes the following tools:
+ KMC2 for k-mer counting,
+ kmc_tools for trimming a k-mer database.

### Compilation
RECKONER is available for Linux. To compile RECKONER
use g++ version 4.9.2 or newer by typing "make",
executable files will be created in the "bin" directory.

### Running
To run RECKONER execute script "run.sh" from the "bin" directory (typically
by typing "./run.sh ARGUMENTS").

As ARGUMENTS type the following:
+ k (k-mer length),
+ output directory; if it does not exist, it will be created,
+ number of threads,
+ list of input FASTQ files.

After running, RECKONER will perform the following actions:
+ k-mer counting with KMC2,
+ determining of a cutoff threshold,
+ k-mer database trimming and converting,
+ read error correcting.

Result files will be created in the output directory with a "corrected" suffix.

Examples:
```
./run.sh 25 results 16 SRR088579.fastq.gz
./run.sh 30 results 8 ERR729973.fastq ERR729974.fastq
```

### Miscellaneous
+ Homepage: http://sun.aei.polsl.pl/REFRESH/
+ Contact: maciej.dlugosz@polsl.pl
+ RECKONER is distributed under license GPL 3.
