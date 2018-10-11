DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
mkdir sortbyname -p
xargs -a bamlist.txt -n1 -P 10 -I FILE $DIR/view_xargs.sh FILE 
