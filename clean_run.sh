#!/bin/bash
if [ "$1" == "" ]
then
	echo Please enter a directoryb to clean..
	echo "Command format clean_run.sh <directory name>"
	exit 1
fi

rm -r  $1/bed
rm -r $1/last
rm -r $1/clusterlist.opt.subopt.txt
rm -r $1/cluster.opt.subopt
rm -r $1/combined
rm -r $1/contigsx
rm -r $1/convert2fastq.log
rm -r $1/fasta
rm -r $1/blast
rm -r $1/fastq
rm -r $1/final_hits.txt
rm -r $1/HBVALL
rm $1/head.fa.blast
rm $1/head.hbv.txt
rm $1/head.hg19.opt.subopt.txt
rm $1/head.mis.opt.subopt.txt
rm -r $1/log
rm -r $1/oneside
rm -r $1/out
rm -r $1/predictions.msa
rm $1/predictions.opt.subopt.txt*
rm -r $1/readsx
rm -r $1/sortbyname
rm $1/tail.fa.blast
rm $1/tail.hbv.txt
rm $1/tail.hg19.opt.subopt.txt
rm $1/tail.mis.opt.subopt.txt
rm -r $1/tmp.batvi
rm $1/t.opt.subopt.cluster
rm $1/t.opt.subopt.sort
rm $1/t.uniq.clusterx
rm -r $1/unbug
rm -r $1/unmapped
rm -r $1/batvirun.log
rm -r $1/run.log
rm -r $1/*.rc.txt*
