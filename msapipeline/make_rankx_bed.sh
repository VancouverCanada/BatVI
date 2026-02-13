set -o pipefail
set -e
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
Target=$1
Path=${Target##*/}
mkdir $1/bed -p
tr : '\t' < $1/t.opt.subopt.cluster |awk '{St=$2-1;Ed=$2+1;print $1"\t"St"\t"Ed}'|sort -u > $1/bed/$Path.rankx.bed ||exit 1
