DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
echo Processing $1
$DIR/unbug.pl sortbyname/$1.sam unbug/$1.sam $2 $3 2>sortbyname/$1.unpaired
