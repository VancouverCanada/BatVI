#/bin/bash
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
while read F
do
echo Processing $F...
$DIR/msapipeline/call_bpx.sh $F
#$DIR/msapipeline/collate_results.pl $F tmp.batvi/$F.predictionsx.msa.txt $F/predictions.txt >$F.collate.txt 
done<$1
