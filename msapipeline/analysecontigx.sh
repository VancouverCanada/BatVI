set -o pipefail
set -e
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
Target=$1
Path=${Target##*/}
#$DIR/cluster.pl $1/readsx/*.fa|$DIR/merge_msa.pl  >$1/$Path.rc.txt
#fix for large files..
$DIR/cluster.pl $1/readsx |$DIR/merge_msa.pl  >$1/$Path.rc.txt
java -cp $DIR msa.Main $1/$Path.rc.txt 
echo $1/$Path.rc.txt 
