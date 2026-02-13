# BatVI Progress Checklist

Last updated: 2026-02-13

## Current status

- [x] Repository imported and version-controlled.
- [x] Core shell/perl scripts present in root directory.
- [x] Added preflight validator: `preflight_check.sh`.
- [x] Fixed main entry script argument handling: `call_integrations.sh`.
- [x] Fixed build script to fail fast on missing source tree: `build.sh`.
- [x] Corrected README command drift (`gen_paths.sh`, `call_integrations.sh`).
- [x] Restored missing source directories from `brownmp/batvi:devel` image layer (`sha256:e010d62dbdc504076bef816e056cc4be838702954ac1c6232fd9c1d6c6d27112`).
- [x] Added Ubuntu Docker packaging (`docker/Dockerfile`, `docker/entrypoint.sh`).
- [x] Added Picard compatibility wrapper for modern environments (`picard_samtofastq.sh`).
- [x] Added migration-facing docs to Veclon (`README.md`).
- [x] Added bilingual deployment documentation (`README.md`, `docker/README.md`).
- [x] Added Docker-friendly config template (`batviconfig.template.txt`).

## Blocking items (must complete before pipeline can run)

- [x] Restore missing source directories:
  - `BatMis-3.00`
  - `batindel`
  - `bin`
  - `msapipeline`
  - `test` (directory restored; sample data files are not present)
- [x] Install runtime dependencies (inside Docker image):
  - `blastn`
  - `bwa`
  - `samtools`
  - `bedtools`
  - Java runtime (`java -version` must work)
- [ ] Create/verify `batviconfig.txt`.
- [ ] Prepare human/pathogen indexes referenced by config.
- [ ] Prepare processing input directory with `filelist.txt` and FASTQ files.

## Recommended run sequence

1. `bash build.sh`
2. `bash preflight_check.sh <processing_directory>`
3. `bash call_integrations.sh <processing_directory> -t <threads>`

If preflight fails, fix each `[FAIL]` item first and rerun.
