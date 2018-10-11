set -e
set -u
set -o pipefail

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd $1
if [ -e "log/step2" ]
then 
	echo "Step 2 seems complete. Moving to next step..";cd $OLDPWD 
else
	mkdir fastq -p
	cat out/*_1.fa >fastq/head.fa
	cat out/*_2.fa >fastq/tail.fa

	T1=$((Threads/2))
	if (( T1 == 0 ))
	then
		T1=1
	fi
	T2=$((Threads-T1))
	echo -e "Threading: $T1 and $T2\n"

	mkdir blast -p
	${BLAST_PATH}/blastn -db $PATHOGEN_BLAST_DB -query fastq/head.fa -num_threads $T1 -num_descriptions 1 -num_alignments 1 |$DIR/bin/blast_parser.pl > blast/head.hbv.txt &  
	if (( T2 != 0 ))
	then
		${BLAST_PATH}/blastn -db $PATHOGEN_BLAST_DB -query fastq/tail.fa -num_threads $T2 -num_descriptions 1 -num_alignments 1 |$DIR/bin/blast_parser.pl > blast/tail.hbv.txt &  
		wait $!; Err=$?;
    if [ $Err -ne 0 ]; then echo "Blast Error $Err";exit 1;fi 
	else
		wait $!; Err=$?;
    if [ $Err -ne 0 ]; then echo "Blast Error $Err";exit 1;fi 
		${BLAST_PATH}/blastn -db $PATHOGEN_BLAST_DB -query fastq/tail.fa -num_threads 1 -num_descriptions 1 -num_alignments 1 |$DIR/bin/blast_parser.pl > blast/tail.hbv.txt  
    if [ $? -ne 0 ]; then echo "Blast Error $?";exit 1;fi 
	fi
	touch log/step2
	cd $OLDPWD 
fi
