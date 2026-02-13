set -o pipefail
set -e
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
Target=$1
Path=${Target##*/}
HOST=`hostname`

mkdir $1/tmp.batvi -p
$DIR/make_rankx_bed.sh $1 || { echo make_rankx_bed.sh failed;exit 1; }
. $DIR/rankx_hits.sh $1 || { echo rankx_hits.sh failed..;exit 1;}
$DIR/analysecontigx.sh $1 || { echo analysecontigx.sh failed..;exit 1; }
mkdir $1/contigsx/blast -p
mkdir $1/fasta -p
mkdir $1/predictions.msa -p
$DIR/parse_msa_out.pl $1/$Path.rc.txt.out|$DIR/filter_rc.pl >$1/tmp.batvi/$Path.x.fa || { echo parse_msa_out.pl failed;exit 1; }

${BWA_PATH}/bwa aln -t $Threads $PATHOGEN_BWA $1/tmp.batvi/$Path.x.fa >$1/tmp.batvi/$Path.sa || { echo BWA alignment to $PATHOGEN_BWA failed..;exit 1; }
${BWA_PATH}/bwa samse $PATHOGEN_BWA $1/tmp.batvi/$Path.sa $1/tmp.batvi/$Path.x.fa >$1/tmp.batvi/$Path.hbv.sam || { echo BWA to $PATHOGEN_BWA failed..;exit 1; }
${BWA_PATH}/bwa aln -t $Threads $HG_BWA $1/tmp.batvi/$Path.x.fa >$1/tmp.batvi/$Path.sa || { echo BWA alignment to $HG_BWA failed..;exit 1; }
${BWA_PATH}/bwa samse $HG_BWA $1/tmp.batvi/$Path.sa $1/tmp.batvi/$Path.x.fa >$1/tmp.batvi/$Path.hg.sam || { echo BWA alignment to $HG_BWA failed..;exit 1; } 
cat $1/tmp.batvi/$Path.*.sam |sort -k1,1 >$1/tmp.batvi/$Path.allmaps.sam
$DIR/remove_full_maps.pl $1/tmp.batvi/$Path.allmaps.sam >$1/contigsx/$Path.fa || { echo remove_full_maps.pl failed..;exit 1; }

time ${BLAST_PATH}/blastn -db $HG_BLAST_DB -query $1/contigsx/$Path.fa -num_threads $Threads >$1/contigsx/blast/$Path.hg.blast || { echo blast to $HG_BLAST_DB failed..;exit 1; }

$DIR/blast_parserX.pl $1/contigsx/blast/$Path.hg.blast >$1/contigsx/blast/$Path.hg.txt || { echo Error parsing blast hits..;exit 1; }
$DIR/get_probablyhbv.pl $1/contigsx/$Path.fa $1/contigsx/blast/$Path.hg.txt >$1/fasta/$Path.splitmap.fa || { echo Error extracting potential HBV hits..;exit 1; }
awk '(NR%2==0)' $1/fasta/$Path.splitmap.fa |sort -u |awk '{print ">"$0"\n"$0}'>$1/fasta/$Path.minimal.fa || { echo Error Generating fastat $Path.minimal.fa;exit 1; }
${BLAST_PATH}/blastn -db $PATHOGEN_BLAST_DB -query $1/fasta/$Path.minimal.fa -num_threads $Threads -word_size 20 |$DIR/blast_parser.pl >$1/contigsx/blast/$Path.minimal.vir.txt || { echo Error blasting to $PATHOGEN_BLAST_DB;exit 1; }
$DIR/print_bpx.pl $1/contigsx/blast/$Path.minimal.vir.txt $1/fasta/$Path.splitmap.fa >$1/tmp.batvi/$Path.predictionsx.msa.txt || { echo Cannot report breakpoints;exit 1; }

