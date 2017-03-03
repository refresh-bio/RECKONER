## RECKONER - Read Error Corrector Based on KMC
### RECKONER
RECKONER is a tool for correction of Illumina reads. It is based
on the error correction algorithm BLESS and utilizes the following tools:
+ KMC2 for k-mer counting,
+ KMC tools for trimming a k-mer database.

RECKONER is available for Linux and Windows.

### Compilation
To compile RECKONER in Linux use g++ version 4.9.2 or newer
by typing "make" in the project directory, as a result
there will be generated executable files in the "bin" directory.

To compile RECKONER in Windows open file "RECKONER.sln"
in Visual Studio 2013, select "Release" from the drop-down list.
Then click "BUILD" and "Build Solution" from the main menu.

### Running
To run RECKONER in Linux execute file "reckoner" from the "bin" directory 
(typically by typing "./reckoner ARGUMENTS").

To run RECKONER in Windows execute file "reckoner.exe" from the "bin" directory
(typically by typing "reckoner.exe ARGUMENTS" from the command line).

As ARGUMENTS type the following:
+ -help - prints short help,
+ -read FASTQ_FILE - input FASTQ read file name (possibly gzipped), can be passed many times,
+ -prefix DIRECTORY - directory for temporary and output files, (optional, default current directory (.)),
+ -kmerlength K - length of k-mers (optional),
+ -genome G - approximate genome size (optional),
+ -memory N - max k-mer counting memory consumption in GB (optional, default 4),
+ -extend N - max extend length, (optional, default 2),
+ -threads N - number of correcting threads, (optional, default number of available virtual cores).

If k-mer length is not given it is determined automatically.
Passing the approximate genome size may slightly improve
k-mer length calculation.

Result files will be created in the output directory with a "corrected" suffix.

Instead of using -read option the input FASTQ files
can be specified at the end of the command line.

Examples:
```
./reckoner -kmerlength 25 -prefix results -threads 16 -read SRR088579.fastq.gz
./reckoner -read ERR729973.fastq -read ERR729974.fastq
./reckoner -kmerlength 30 ERR729973.fastq ERR729974.fastq
```

### Miscellaneous
+ Homepage: http://sun.aei.polsl.pl/REFRESH/
+ Contact: maciej.dlugosz@polsl.pl
+ RECKONER is distributed under license GPL 3.
