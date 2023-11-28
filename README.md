## RECKONER - Read Error Corrector Based on KMC
### RECKONER
RECKONER is a tool for Illumina reads correction. It is based
on the algorithm BLESS, supplies indels correction and utilizes the following tools:
+ KMC3 for k-mer counting,
+ KMC tools for trimming a k-mer database.

RECKONER is available for Linux and Windows.

### Building
Building RECKONER in Linux requires presence of Python 3 libraries
(install python3-dev package in Debian or python3-devel package in Fedora).
To compile use g++ 6.3.0 or newer by typing "make -j8"
in the project directory. As a result
executable files will be generated in a "bin" directory.

To compile RECKONER in Windows open file "RECKONER.sln"
in Visual Studio 2022 or newer, select "Release" from the drop-down list.
Then click "Build" and "Build Solution" from the main menu.

You can also use precompiled binaries. In Windows it could
be needed to download Visual C++ 2015 Redistributable.

### Running
To run RECKONER in Linux execute a file "reckoner" from the "bin" directory 
(typically by typing "./reckoner {ARGUMENTS} {FASTQ files}") from the command line.

To run RECKONER in Windows execute a file "reckoner.exe" from the "bin" directory
(typically by typing "reckoner.exe {ARGUMENTS} {FASTQ files}") from the command line.

In both the cases the supporting tools (i.e. kmc(.exe) and kmc_tools(.exe)) have to be
present in the same directory.

A list of allowed {ARGUMENTS} is given below.
+ Main arguments (sufficient in a general use):
	+ -help - print a short help,
	+ -read <FASTQ file> - input FASTQ(.gz) read file name, can be passed many times,
	+ -kmerlength <K> - k-mers length, not mandatory,
+ Supporting:
	+ -genome <G> - approximate genome size, not mandatory - when -kmerlength is not given, it can ease k-mer length determination,
+ Quality tunning:
	+ -longkmer - verify corrections with long k-mers (worthful for _de novo_ assembly, requires more RAM and CPU time), it enables:
		+ -accept - if it is not possible to verify some read, keep the best of its possible corrections (default: do not correct the read),
		+ -longratio <RATIO> - long to normal k-mer length ratio, must be > 1.0,
	+ -noindels - do not correct insertions and deletions (substitutions correction only),
	+ -extend N - max extension length,
+ Formal issues:
	+ -prefix <DIRECTORY> - output directory (optional, by default the current directory (.)),
	+ -threads <N> - number of correcting threads (optional, by default number of available virtual cores),
	+ -mark - mark corrected reads by adding an information to their headers,
	+ -fixlength - if the read headers contain read length (e.g. \"length=100\"), fix this value after indel correction (cannot be used with -noindels),
+ Useful for debugging:
	+ -kmcmemory <N> - k-mer counting memory consumption limit in GB,
	+ -kmcram - count k-mers in RAM-only mode,
	+ -reuse - do not remove a KMC database; if a proper database exists, use it,
	+ -verbose - return more statistics.


Instead of using -read option the input FASTQ files
can be specified at the end of the command line.

K-mer length, if not given, will be determined automatically.
Delivering genome size can facilitate and accelerate this determination.

To correct paired-end files simply pass both of them
as two single-end files. RECKONER will keep the pairs in files.

Result files will be created in the output directory with a "corrected" suffix
in the <DIRECTORY> directory.

Examples:
```
./reckoner -kmerlength 25 -prefix results -threads 16 -read SRR088579.fastq.gz
./reckoner -read ERR729973.fastq -read ERR729974.fastq
./reckoner -kmerlength 30 ERR729973.fastq ERR729974.fastq
./reckoner -kmerlength 20 ERR422544_1.fastq.gz ERR422544_2.fastq.gz
```
To achieve reads best for _de novo_ assembly:
```
./reckoner -kmerlength 20 -longkmer ERR422544_1.fastq.gz ERR422544_2.fastq.gz
```

### Debugging in Windows
To debug RECKONER in Windows, select "Debug" from the Visual Studio
drop-down list. In the RECKONER project properties set
the working directory as:
```
$(SolutionDir)bin\$(Configuration)
```
Input files should be located in "bin\Debug".

### Miscellaneous
+ Contact: maciej.dlugosz@polsl.pl
+ RECKONER is distributed under license GPL 3.

## Citing
+ <a href="https://doi.org/10.1007/978-3-319-67792-7_12">Długosz, M., Deorowicz, S., Kokot, M. (2017) Improvements in DNA Reads Correction, International Conference on Man-Machine Interactions, 115&ndash;124.</a>
+ <a href="https://doi.org/10.1093/bioinformatics/btw746">Długosz, M., Deorowicz, S. (2017) RECKONER: read error corrector based on KMC, Bioinformatics, 33:1086&ndash;1089.</a>