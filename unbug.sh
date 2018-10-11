DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
mkdir unbug -p
rm unbug/*
xargs -a bamlist.txt -n1 -P 2 $DIR/unbug_xarg.sh 
xargs -a bamlist.txt -n1 -P 2 $DIR/unbug_xarg.sh 
#../relabelunbugreads.sh 
