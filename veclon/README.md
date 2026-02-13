# AccenTrust Veclon

Veclon is AccenTrust's secure orchestration platform for viral integration analysis across heterogeneous engines.

## Platform modules

- `Veclon Core`: API, workflow orchestration, policy and audit.
- `Veclon Bridge`: external engine adapters (user-supplied tools only).
- `Veclon Studio`: visualization and operations UI.
- `Veclon Copilot`: local AI assistant for interpretation and reporting.

## Closed-source-friendly boundary

- Veclon does not bundle GPL engines or source code.
- Engines are mounted or referenced externally by path.
- Bridge adapters execute external commands and normalize outputs into canonical schemas.

## Repository map

- `apps/`: deployable services (`core-api`, `worker`, `studio`, `copilot`)
- `packages/`: shared libraries, adapters, fusion and reporting modules
- `contracts/`: JSON schemas for engine manifests and canonical outputs
- `infra/`: docker/compose/k8s deployment assets
- `docs/`: architecture, API, compliance docs

## Quick start (scaffold stage)

1. Review `contracts/*.schema.json`.
2. Start infra dependencies with `infra/compose/docker-compose.yml`.
3. Run Core API skeleton from `apps/core-api`.
