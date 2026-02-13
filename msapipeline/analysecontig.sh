set -e
set -u
set -o pipefail

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
$DIR/cluster.pl $1/reads/*.fa >$1.rc.txt || exit 1
java -cp $DIR msa.Main $1.rc.txt  || exit 1
