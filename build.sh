#!/usr/bin/env bash
set -euo pipefail

CUR_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

usage() {
  cat <<'EOF'
Usage:
  bash build.sh [clean]
EOF
}

if [[ $# -gt 1 ]]; then
  usage
  exit 2
fi

if [[ "${1:-}" == "clean" ]]; then
  [[ -d "$CUR_DIR/BatMis-3.00" ]] && (cd "$CUR_DIR/BatMis-3.00" && make clean)
  [[ -d "$CUR_DIR/batindel" ]] && (cd "$CUR_DIR/batindel" && make clean)
  if [[ -d "$CUR_DIR/bin" ]]; then
    rm -f "$CUR_DIR/bin/cluster_bp" "$CUR_DIR"/bin/*.o
  fi
  if [[ -d "$CUR_DIR/msapipeline/msa" ]]; then
    rm -f "$CUR_DIR"/msapipeline/msa/*.class
  fi
  exit 0
fi

for d in BatMis-3.00 batindel bin msapipeline/msa; do
  if [[ ! -d "$CUR_DIR/$d" ]]; then
    echo "Missing required directory: $CUR_DIR/$d" >&2
    echo "Please restore the full BatVI source tree before building." >&2
    exit 1
  fi
done

echo "=========================================================="
echo "               COMPILING BatMis"
echo "=========================================================="
(cd "$CUR_DIR/BatMis-3.00" && ./configure && make && make copy)

echo "=========================================================="
echo "               COMPILING BATINDEL-lite"
echo "=========================================================="
(cd "$CUR_DIR/batindel" && ./configure && make)

echo "=========================================================="
echo "               COMPILING Binaries"
echo "=========================================================="
(cd "$CUR_DIR/bin" && bash ./build.sh)

echo "=========================================================="
echo "               COMPILING MSA"
echo "=========================================================="
(cd "$CUR_DIR/msapipeline/msa" && javac *.java)

bash "$CUR_DIR/manualcompile.sh"

command -v blastn >/dev/null 2>&1 || { echo >&2 "Please check if BLAST is installed"; }
command -v bwa >/dev/null 2>&1 || { echo >&2 "Please check if BWA is installed"; }
