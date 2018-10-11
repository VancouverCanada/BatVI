set -e
set -u
set -o pipefail

CUR_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

>&2 echo -e "[\033[0;31m  Loading configuration information \033[0m ]"
if [ -e $1/batviconfig.txt ] 
then
	. $1/batviconfig.txt
else
	if [ ! -e "$CUR_DIR/batviconfig.txt" ]; then >&2 echo "Configuration file is missing..";exit 1;fi
	. $CUR_DIR/batviconfig.txt
fi

>&2 echo -e "[\033[0;31m  Extracting potential Pathogen reads \033[0m ]"
. $CUR_DIR/commandline.sh "$@" ||return 1
. $CUR_DIR/extract_hbv_from_fasta.sh $1 ||return 1
####$CUR_DIR/get_blast_hits.sh $1
>&2 echo -e "[\033[0;31m  Searching for reads with pathogen presence \033[0m ]"
. $CUR_DIR/hbvblast.sh $1 || return 1
>&2 echo -e "[\033[0;31m  Calling integrations \033[0m ]"
. $CUR_DIR/bin/integrations.sh $1 ||return 1
>&2 echo -e "[\033[0;31m  Local assembling clusters\033[0m ]"
. $CUR_DIR/msapipeline/call_bpx.sh $1 ||return 1
>&2 echo -e "[\033[0;31m  Report rank 1 integrations \033[0m ]"
$CUR_DIR/combine_hits.pl $1 200 > $1/final_hits.txt

