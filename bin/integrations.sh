set -e
set -u
set -o pipefail

trap killgroup SIGINT

killgroup()
{
	  echo killing...
	  kill 0
}

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd $1
LIB=$1
HOST=`hostname`

T1=$((Threads/2))
if (( T1 == 0 ))
then
	T1=1
fi

T2=$((Threads-T1))
echo -e "Threading: $T1 and $T2\n"


if [ ! -e "log/cluster" ]
then
	mkdir HBVALL/fastq -p
	mkdir HBVALL/tmp -p
	ln -s $PWD/fastq/head.fa HBVALL/fastq/unmapped_1.fa ||true
	ln -s $PWD/fastq/tail.fa HBVALL/fastq/unmapped_2.fa ||true
	ln -s $PWD/blast/*.txt . ||true
	$DIR/get_hbv_list.sh HBVALL || { echo get_hbv_list.sh failed.. ; exit 1; } 
	$DIR/extract_unmapped_pairsX.pl HBVALL || { echo extract_unmapped_pairsX.pl failed.. ; exit 1; }
	echo "time ${BLAST_PATH}/blastn -db $HG_BLAST_DB -query HBVALL/tmp/tail.fq -num_threads $T1" 
	time ${BLAST_PATH}/blastn -db $HG_BLAST_DB -query HBVALL/tmp/tail.fq -num_threads $T1 >tail.fa.blast &
	if (( T2 != 0 ))
	then
    if [ $? -ne 0 ]; then echo "Blast Error $?";exit 1;fi 
		time ${BLAST_PATH}/blastn -db $HG_BLAST_DB -query HBVALL/tmp/head.fq -num_threads $T2 >head.fa.blast &
		wait $!; Err=$?;
    if [ $Err -ne 0 ]; then echo "Blast Error $Err";exit 1;fi 
	else
		wait $!; Err=$?;
    if [ $Err -ne 0 ]; then echo "Blast Error $Err";exit 1;fi 
		time ${BLAST_PATH}/blastn -db $HG_BLAST_DB -query HBVALL/tmp/head.fq  >head.fa.blast 
	fi
  if [ $? -ne 0 ]; then echo "Blast Error $?";exit 1;fi 

	$DIR/blast_parserX.pl tail.fa.blast|$DIR/pick_opt_subopt.pl >tail.hg19.opt.subopt.txt || { echo blast_parserX.pl tail.fa.blast failed.. ; exit 1; } 
	$DIR/blast_parserX.pl head.fa.blast|$DIR/pick_opt_subopt.pl >head.hg19.opt.subopt.txt || { echo blast_parserX.pl head.fa.blast failed.. ; exit 1; }


	$DIR/parse_blast_2_clust.pl head.hg19.opt.subopt.txt tail.hg19.opt.subopt.txt |sort -k1,1 -k2,2n >t.opt.subopt.sort || { echo parse_blast_2_clust.pl failed.. ; exit 1; } 
	$DIR/cluster.pl t.opt.subopt.sort >t.opt.subopt.cluster || { echo cluster.pl failed.. ; exit 1; }
	touch log/cluster
fi
awk -F: '{print $1"\t"$2}' t.opt.subopt.cluster |sort -u >clusterlist.opt.subopt.txt || { echo cluster list collection failed.. ; exit 1; }

mkdir -p cluster.opt.subopt 
$DIR/split.pl clusterlist.opt.subopt.txt cluster.opt.subopt/cluster $Threads || { echo split.pl failed.. ; exit 1; } 
echo Finding Breakpoints..
for F in cluster.opt.subopt/*
do
	N=${F##*.}
	while read L;
	do 
	$DIR/cluster_bp $L head.hg19.opt.subopt.txt tail.hg19.opt.subopt.txt head.hbv.txt tail.hbv.txt $HG_GENOME.pac HBVALL/tmp/head.fq HBVALL/tmp/tail.fq ${HG_GENOME}.ann.location $Log_File.$N $Dup_Filter 
	done<$F > predictions.opt.subopt.txt.$N &
done
wait $!; Err=$?;
if [ $Err -ne 0 ]; then echo "Blast Error $Err";exit 1;fi 

echo Collecting hits..
cat predictions.opt.subopt.txt.* >predictions.opt.subopt.txt
if [ "$Log_File" != "" ]
then
	cat $Log_File.* >$Log_File || return 1
	rm $Log_File.* 
fi

sort -k8,8n -k2,2 -k3,3n -u predictions.opt.subopt.txt >xxx || { echo sorting predictions failed.. ; exit 1; } ;mv xxx predictions.opt.subopt.txt
cd -
