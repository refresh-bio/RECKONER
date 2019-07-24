## RECKONER - Read Error Corrector Based on KMC
### RECKONER
RECKONER is a tool for correction of Illumina reads. It is based
on the error correction algorithm BLESS and utilizes the following tools:
+ KMC3 for k-mer counting,
+ KMC tools for trimming a k-mer database.

RECKONER is available for Linux and Windows.

### Building
Building RECKONER in Linux requires presence of Python 3 libraries
(install python3-dev package in Debian or python3-devel package in Fedora).
To compile use g++ 4.9.2 or newer
by typing "make" in the project directory. As a result
executable files will be generated in the "bin" directory.

To compile RECKONER in Windows open file "RECKONER.sln"
in Visual Studio 2015, select "Release" from the drop-down list.
Then click "Build" and "Build Solution" from the main menu.

You can also use precompiled binaries. In Windows there could
be need to download Visual C++ 2015 Redistributable.

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

Instead of using -read option the input FASTQ files
can be specified at the end of the command line.

To correct paired-end files simply pass both of them
as two single-end files. RECKONER will keep the pairs in files.

Result files will be created in the output directory with a "corrected" suffix.

Examples:
```
./reckoner -kmerlength 25 -prefix results -threads 16 -read SRR088579.fastq.gz
./reckoner -read ERR729973.fastq -read ERR729974.fastq
./reckoner -kmerlength 30 ERR729973.fastq ERR729974.fastq
./reckoner -kmerlength 20 ERR422544_1.fastq.gz ERR422544_2.fastq.gz
```

### Miscellaneous
+ Homepage: http://sun.aei.polsl.pl/REFRESH/
+ Contact: maciej.dlugosz@polsl.pl
+ RECKONER is distributed under license GPL 3.
