#!/bin/bash
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd $1

HOST=`hostname`
LIB=$1
if [ "$2" == "" ]
then
	if [ ! -e "bamlist.txt" ]
	then
		ls fastq/*_1.fastq.gz|sed 's/fastq\///' >bamlist.txt
	fi
	rm fastq/*
	. $DIR/extract_unmapped_and_oneside.sh
	xargs -a bamlist.txt -n1 -P 2 -I FILE $DIR/extracthbv.sh sam/FILE

	. $DIR/join_sam.sh
	. $DIR/sort_by_name.sh
	. $DIR/cat_sorted_sam.sh
fi
. $DIR/unbug.sh
. $DIR/convert_to_fastq.sh unbug &>convert2fastq.log



