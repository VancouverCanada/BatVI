#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

printf "[Veclon] root: %s\n" "$ROOT_DIR"
printf "[Veclon] starting local infra...\n"

docker compose -f "$ROOT_DIR/infra/compose/docker-compose.yml" up -d

printf "[Veclon] infra started.\n"
printf "- Postgres: localhost:5432\n"
printf "- Redis: localhost:6379\n"
printf "- MinIO API: localhost:9000\n"
printf "- MinIO Console: localhost:9001\n"
