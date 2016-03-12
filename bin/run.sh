#!/bin/bash

#arguments:
#1 - k
#2 - temp and output directory
#3 - number of threads
#4... - input files list

if [ $# -lt 4 ]
then
	echo "Error: not enough parameters given." >&2
	echo "Try:" >&2
	echo "./run.sh k tmp_dir threads FASTQ_file1 FASTQ_file2 ..." >&2
	exit
fi

trap 'rm -fr $4.lst $4.lst.tmp.kmc_pre $4.lst.tmp.kmc_suf $4.lst.kmc_pre $4.lst.kmc_suf' EXIT

mkdir -p $2

echo $4 > $4.lst
file_list="-read $4"

for ((i=5; $i <= $#; i++))
do
	echo ${!i} >> $4.lst
	file_list="$file_list -read ${!i}"
done

echo '##################################################'
echo K-MER COUNTING
echo '##################################################'
./kmc -k$1 -m4 -ci2 @$4.lst $4.lst.tmp $2
rm $4.lst
if [ ! -e $4.lst.tmp.kmc_pre ]
then
	echo "Error: failed to execute KMC." >&2
	exit
fi
echo DONE
echo

echo '##################################################'
echo DETERMINING CUTOFF THRESHOLD
echo '##################################################'
cut_res=`./cutter $4.lst.tmp`
if [ ! $? -eq 0 ]
then
	echo "Error: failed to determine the cutoff threshold." >&2
	exit
fi
echo "Cutoff: $cut_res"
echo DONE
echo

echo '##################################################'
echo REMOVING UNTRUSTED K-MERS
echo '##################################################'
./kmc_tools reduce $4.lst.tmp -ci$cut_res $4.lst
if [ ! $? -eq 0 ]
then
	echo "Error: failed to remove distrusted k-mers." >&2
	exit
fi
echo DONE
echo

echo '##################################################'
echo CORRECTING ERRORS
echo '##################################################'
./reckoner -kmerlength $1 $file_list -prefix $2 -threads $3
if [ ! $? -eq 0 ]
then
	echo "Error: failed to correct errors." >&2
	exit
fi
echo DONE
echo
rm $4.lst.kmc_pre $4.lst.kmc_suf

