#!/usr/bin/env bash
set -euo pipefail

CUR_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

usage() {
  cat <<'EOF'
Usage:
  bash preflight_check.sh <processing_directory>

Checks whether BatVI has the files, directories, tools, and config needed
before running call_integrations.sh.
EOF
}

if [[ $# -ne 1 ]]; then
  usage
  exit 2
fi

PROCESSING_DIR="$1"
if [[ ! -d "$PROCESSING_DIR" ]]; then
  echo "[FAIL] Processing directory not found: $PROCESSING_DIR"
  exit 1
fi

FAILS=0
WARNS=0

pass() { echo "[PASS] $1"; }
warn() { echo "[WARN] $1"; WARNS=$((WARNS+1)); }
fail() { echo "[FAIL] $1"; FAILS=$((FAILS+1)); }

require_dir() {
  local d="$1"
  if [[ -d "$d" ]]; then
    pass "Directory exists: $d"
  else
    fail "Missing directory: $d"
  fi
}

require_file() {
  local f="$1"
  if [[ -e "$f" ]]; then
    pass "File exists: $f"
  else
    fail "Missing file: $f"
  fi
}

require_cmd() {
  local c="$1"
  if command -v "$c" >/dev/null 2>&1; then
    pass "Command available: $c"
  else
    fail "Command missing: $c"
  fi
}

check_java_runtime() {
  if command -v java >/dev/null 2>&1 && java -version >/dev/null 2>&1; then
    pass "Java runtime is available"
  else
    fail "Java runtime is not available"
  fi
}

check_non_empty_var() {
  local var_name="$1"
  local v="${!var_name:-}"
  if [[ -n "$v" ]]; then
    pass "Config value set: $var_name"
  else
    fail "Config value missing: $var_name"
  fi
}

check_bwa_index() {
  local label="$1"
  local base="$2"
  local ok=1
  for ext in amb ann bwt pac sa; do
    if [[ ! -e "${base}.${ext}" ]]; then
      ok=0
      break
    fi
  done
  if [[ $ok -eq 1 ]]; then
    pass "BWA index files present for $label: $base.[amb|ann|bwt|pac|sa]"
  else
    fail "BWA index files incomplete for $label: $base"
  fi
}

check_batmis_index() {
  local base="$1"
  if [[ -e "${base}.pac" ]] && ([[ -e "${base}.ann.location" ]] || [[ -e "${base}.ann.locations" ]]); then
    pass "BatMis index files present: ${base}.pac and .ann.location(s)"
  else
    fail "BatMis index files incomplete for INDEX: $base"
  fi
}

check_blast_db() {
  local label="$1"
  local base="$2"
  if [[ -e "${base}.nhr" ]] && [[ -e "${base}.nin" ]] && [[ -e "${base}.nsq" ]]; then
    pass "BLAST db present for $label: $base.[nhr|nin|nsq]"
  else
    fail "BLAST db files incomplete for $label: $base"
  fi
}

check_bin_path() {
  local label="$1"
  local dir="$2"
  local cmd="$3"
  if [[ -x "${dir}/${cmd}" ]]; then
    pass "Binary path valid: ${label} (${dir}/${cmd})"
  else
    fail "Binary path invalid for ${label}: expected executable ${dir}/${cmd}"
  fi
}

echo "== BatVI Preflight =="
echo "Workspace: $CUR_DIR"
echo "Processing dir: $PROCESSING_DIR"

echo
echo "## Repository layout"
require_dir "$CUR_DIR/BatMis-3.00"
require_dir "$CUR_DIR/batindel"
require_dir "$CUR_DIR/bin"
require_dir "$CUR_DIR/msapipeline"
require_file "$CUR_DIR/call_integrations.sh"
require_file "$CUR_DIR/combine_hits.pl"

echo
echo "## Runtime commands"
require_cmd blastn
require_cmd bwa
require_cmd samtools
require_cmd bedtools
require_cmd perl
check_java_runtime

echo
echo "## Processing directory"
require_file "$PROCESSING_DIR/filelist.txt"

CFG_FILE="$PROCESSING_DIR/batviconfig.txt"
if [[ ! -e "$CFG_FILE" ]]; then
  CFG_FILE="$CUR_DIR/batviconfig.txt"
fi
if [[ ! -e "$CFG_FILE" ]]; then
  fail "Configuration file not found (checked $PROCESSING_DIR/batviconfig.txt and $CUR_DIR/batviconfig.txt)"
else
  pass "Using config: $CFG_FILE"

  # shellcheck disable=SC1090
  . "$CFG_FILE"

  echo
  echo "## Config values"
  check_non_empty_var INDEX
  check_non_empty_var PATHOGEN_BLAST_DB
  check_non_empty_var HG_BLAST_DB
  check_non_empty_var HG_GENOME
  check_non_empty_var HG_BWA
  check_non_empty_var PATHOGEN_BWA
  check_non_empty_var BLAST_PATH
  check_non_empty_var BWA_PATH
  check_non_empty_var SAMTOOLS_PATH
  check_non_empty_var BEDTOOLS_PATH

  if [[ -n "${INDEX:-}" ]]; then
    check_batmis_index "$INDEX"
  fi
  if [[ -n "${PATHOGEN_BLAST_DB:-}" ]]; then
    check_blast_db "PATHOGEN_BLAST_DB" "$PATHOGEN_BLAST_DB"
  fi
  if [[ -n "${HG_BLAST_DB:-}" ]]; then
    check_blast_db "HG_BLAST_DB" "$HG_BLAST_DB"
  fi
  if [[ -n "${HG_GENOME:-}" ]]; then
    require_file "$HG_GENOME"
  fi
  if [[ -n "${HG_BWA:-}" ]]; then
    check_bwa_index "HG_BWA" "$HG_BWA"
  fi
  if [[ -n "${PATHOGEN_BWA:-}" ]]; then
    check_bwa_index "PATHOGEN_BWA" "$PATHOGEN_BWA"
  fi

  if [[ -n "${BLAST_PATH:-}" ]]; then
    check_bin_path "BLAST_PATH" "$BLAST_PATH" "blastn"
  fi
  if [[ -n "${BWA_PATH:-}" ]]; then
    check_bin_path "BWA_PATH" "$BWA_PATH" "bwa"
  fi
  if [[ -n "${SAMTOOLS_PATH:-}" ]]; then
    check_bin_path "SAMTOOLS_PATH" "$SAMTOOLS_PATH" "samtools"
  fi
  if [[ -n "${BEDTOOLS_PATH:-}" ]]; then
    check_bin_path "BEDTOOLS_PATH" "$BEDTOOLS_PATH" "bedtools"
  fi
fi

echo
echo "## filelist.txt quick check"
if [[ -e "$PROCESSING_DIR/filelist.txt" ]]; then
  line_no=0
  while IFS= read -r line || [[ -n "$line" ]]; do
    line_no=$((line_no+1))
    [[ -z "$line" ]] && continue
    [[ "$line" =~ ^# ]] && continue
    IFS=';' read -r read1 read2 insertsize <<<"$line"

    if [[ -z "${read1:-}" ]]; then
      fail "filelist.txt line ${line_no}: missing read1"
      continue
    fi

    if [[ -e "$PROCESSING_DIR/$read1" ]] || [[ -e "$read1" ]]; then
      pass "filelist line ${line_no}: read1 found ($read1)"
    else
      fail "filelist line ${line_no}: read1 not found ($read1)"
    fi

    if [[ -n "${read2:-}" ]]; then
      if [[ -e "$PROCESSING_DIR/$read2" ]] || [[ -e "$read2" ]]; then
        pass "filelist line ${line_no}: read2 found ($read2)"
      else
        fail "filelist line ${line_no}: read2 not found ($read2)"
      fi
    else
      warn "filelist line ${line_no}: read2 empty (single-end mode)"
    fi

    if [[ -n "${insertsize:-}" ]] && ! [[ "$insertsize" =~ ^[0-9]+$ ]]; then
      fail "filelist line ${line_no}: insert size is not numeric ($insertsize)"
    fi
  done < "$PROCESSING_DIR/filelist.txt"
fi

echo
echo "== Summary =="
echo "Warnings: $WARNS"
echo "Failures: $FAILS"

if [[ $FAILS -ne 0 ]]; then
  exit 1
fi

echo "Preflight passed."
