set -e
set -u
set -o pipefail
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
mkdir bed -p
awk '($5==1)' $1/predictions.txt |sort -k3,3n|cut -f2,3|$DIR/make_rank1_bed.pl > bed/$1.rank1.bed || return 1
