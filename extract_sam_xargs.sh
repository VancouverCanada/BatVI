#!/bin/bash
# basic pairing..
HOST=`hostname`
EMAIL="drcyber@gmail.com"
EMAILMESSAGE="pair.txt"

HEAD=$1
TAIL=$2
INDEX=$3
LIB=$1
if [ ! -e sam ]
then
	mkdir sam
fi
HEAD=${HEAD#fastq/}
TAIL=${TAIL#fastq/}
${BWA_PATH}/bwa sampe ${INDEX} sai/${HEAD}.sai sai/${TAIL}.sai fastq/${HEAD} fastq/${TAIL} 2>log/${HEAD}.pairing.log |${SAMTOOLS_PATH}/samtools view -bS "-" > sam/${HEAD}.bam 
if [ "$?" != 0 ]
then 
	ZIPBAD=1
	SUBJECT="$HOST:$LIB bwa pairing bad.."
	echo ------------- $2 --------- >> $EMAILMESSAGE
	echo $1 >> $EMAILMESSAGE
	#/bin/mail -s "$SUBJECT" "$EMAIL" < $EMAILMESSAGE
	echo VMAP FAIL  >> log/job.log
	date >>log/job.log
	exit 255
fi
