#/bin/bash
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
F=$1
echo Processing $F...
$DIR/pipeline/call_bpx.sh $F
$DIR/pipeline/collate_results.pl $F tmp/$F.predictionsx.msa.txt $F/predictions.txt >$F.collate.txt 
