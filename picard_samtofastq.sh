#!/usr/bin/env bash
set -euo pipefail

# Backward-compatible SamToFastq launcher:
# 1) legacy SamToFastq.jar under PICARD_PATH
# 2) picard CLI
# 3) distro picard.jar with explicit tool name
if [[ -n "${PICARD_PATH:-}" ]] && [[ -f "${PICARD_PATH}/SamToFastq.jar" ]]; then
  exec java -jar "${PICARD_PATH}/SamToFastq.jar" "$@"
fi

if command -v picard >/dev/null 2>&1; then
  exec picard SamToFastq "$@"
fi

if [[ -f "/usr/share/java/picard.jar" ]]; then
  exec java -jar /usr/share/java/picard.jar SamToFastq "$@"
fi

echo "Cannot find a Picard SamToFastq runtime." >&2
echo "Set PICARD_PATH with SamToFastq.jar, or install picard/picard.jar." >&2
exit 1
