set -o pipefail
set -e
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
Target=$1
Path=${Target##*/}
#Convert cluster list in to bed form, retain only uniq hits..
mkdir $1/tmp.batvi $1/readsx -p
sort -u -k4,4 -k1,1 -k2,2 $1/t.opt.subopt.cluster >$1/t.uniq.clusterx || exit 1
$DIR/cluster2bed.pl $1/t.uniq.clusterx |cut -f1,2,3,6 >$1/tmp.batvi/$Path.t.clusterx || exit 1
#rankx.bed is a bed file consisting of rankx BP, where those BP within 1k are merged
#Extract those reads and intervals corresponding to rankx hits..
${BEDTOOLS_PATH}/bedtools intersect -a $1/bed/$Path.rankx.bed -b $1/tmp.batvi/$Path.t.clusterx -wb|awk '{print $4":"$5":"$6"\t"$7}' > $1/tmp.batvi/$Path.confident.clustersx || exit 1
cp $1/tmp.batvi/$Path.confident.clustersx $1/readsx
$DIR/get_reads_in_clustersx.pl $1/tmp.batvi/$Path.confident.clustersx $1/HBVALL/tmp/head.fq $1/HBVALL/tmp/tail.fq $1 || exit 1
