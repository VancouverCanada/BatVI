DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
mkdir -p "$1/fastq"
xargs -a bamlist.txt -n1 -P 2 -I FILE "$DIR/picard_samtofastq.sh" INPUT="$1/FILE.sam" FASTQ="$1/fastq/FILE_1.fq" SECOND_END_FASTQ="$1/fastq/FILE_2.fq" VALIDATION_STRINGENCY=SILENT
