#!/usr/bin/env bash
set -euo pipefail

CUR_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

usage() {
  cat <<'EOF'
Usage:
  bash call_integrations.sh <processing_directory> [options]

Options:
  -l, --log <log file>
  -t, --threads <thread number>
  -f, --filterdup
  --skip-preflight
EOF
}

if [[ $# -lt 1 ]]; then
  usage
  exit 2
fi

PROCESSING_DIR="$1"
if [[ ! -d "$PROCESSING_DIR" ]]; then
  echo "Cannot find processing directory: $PROCESSING_DIR" >&2
  exit 1
fi

RUN_PREFLIGHT=1
for arg in "$@"; do
  if [[ "$arg" == "--skip-preflight" ]]; then
    RUN_PREFLIGHT=0
  fi
done

if [[ $RUN_PREFLIGHT -eq 1 ]]; then
  bash "$CUR_DIR/preflight_check.sh" "$PROCESSING_DIR"
fi

>&2 echo -e "[\033[0;31m  Loading configuration information \033[0m ]"
if [ -e "$PROCESSING_DIR/batviconfig.txt" ]
then
	. "$PROCESSING_DIR/batviconfig.txt"
else
	if [ ! -e "$CUR_DIR/batviconfig.txt" ]; then >&2 echo "Configuration file is missing..";exit 1;fi
	. $CUR_DIR/batviconfig.txt
fi

>&2 echo -e "[\033[0;31m  Extracting potential Pathogen reads \033[0m ]"
. "$CUR_DIR/commandline.sh" "$@" || exit 1
. "$CUR_DIR/extract_hbv_from_fasta.sh" "$PROCESSING_DIR" || exit 1
####$CUR_DIR/get_blast_hits.sh "$PROCESSING_DIR"
>&2 echo -e "[\033[0;31m  Searching for reads with pathogen presence \033[0m ]"
. "$CUR_DIR/hbvblast.sh" "$PROCESSING_DIR" || exit 1
>&2 echo -e "[\033[0;31m  Calling integrations \033[0m ]"
. "$CUR_DIR/bin/integrations.sh" "$PROCESSING_DIR" || exit 1
>&2 echo -e "[\033[0;31m  Local assembling clusters\033[0m ]"
. "$CUR_DIR/msapipeline/call_bpx.sh" "$PROCESSING_DIR" || exit 1
>&2 echo -e "[\033[0;31m  Report rank 1 integrations \033[0m ]"
"$CUR_DIR/combine_hits.pl" "$PROCESSING_DIR" 200 > "$PROCESSING_DIR/final_hits.txt"
