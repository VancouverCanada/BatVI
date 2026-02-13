DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
F="$1"
F=${F##*/}
echo "processing $F.."
"${SAMTOOLS_PATH}/samtools" view "$1.bam" | grep HBV-B2 > "$F.sam"
"$DIR/picard_samtofastq.sh" INPUT="$F.sam" F="fastq/${F}_1.fq" F2="fastq/${F}_2.fq" VALIDATION_STRINGENCY=SILENT
rm -f "$F.sam"
