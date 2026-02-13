#!/usr/bin/env bash
set -euo pipefail

BATVI_HOME="${BATVI_HOME:-/opt/batvi}"
cd "$BATVI_HOME"

usage() {
  cat <<'EOF'
BatVI container entrypoint

Usage:
  batvi-image shell
  batvi-image build
  batvi-image preflight <processing_directory>
  batvi-image run <processing_directory> [call_integrations options]
  batvi-image <any-command> [args...]
EOF
}

if [[ $# -eq 0 ]]; then
  set -- shell
fi

case "$1" in
  shell)
    shift
    exec /bin/bash "$@"
    ;;
  build)
    shift
    exec bash "$BATVI_HOME/build.sh" "$@"
    ;;
  preflight)
    shift
    if [[ $# -lt 1 ]]; then
      usage
      exit 2
    fi
    exec bash "$BATVI_HOME/preflight_check.sh" "$@"
    ;;
  run)
    shift
    if [[ $# -lt 1 ]]; then
      usage
      exit 2
    fi
    exec bash "$BATVI_HOME/call_integrations.sh" "$@"
    ;;
  -h|--help|help)
    usage
    ;;
  *)
    exec "$@"
    ;;
esac
