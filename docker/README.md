# BatVI Docker Guide / BatVI 容器部署指南

## 1) Purpose / 目的

### English

This folder provides a reproducible Ubuntu runtime for BatVI, intended for macOS users and legacy x86-dependent builds.

### 中文

本目录提供 BatVI 的可复现 Ubuntu 运行环境，重点解决 macOS 与旧版 x86 编译依赖问题。

## 2) Build Image / 构建镜像

```bash
docker build --platform linux/amd64 -f docker/Dockerfile -t batvi:ubuntu .
```

## 3) Container Commands / 容器入口命令

Entrypoint file: `/Users/kenneth/BatVI/docker/entrypoint.sh`

- `shell`: open interactive shell
- `build`: run `build.sh` inside container
- `preflight <dir>`: run `preflight_check.sh`
- `run <dir> [options]`: run `call_integrations.sh`

示例：

```bash
docker run --rm -it --platform linux/amd64 batvi:ubuntu shell
```

## 4) Recommended Data Layout / 推荐数据目录结构

```text
/absolute/path/to/processing
  ├── batviconfig.txt
  ├── filelist.txt
  ├── sampleA_1.fq.gz
  └── sampleA_2.fq.gz
```

`filelist.txt` format:

```text
sampleA_1.fq.gz;sampleA_2.fq.gz;800
```

## 5) Run Preflight / 执行预检

```bash
docker run --rm -it --platform linux/amd64 \
  -v /absolute/path/to/processing:/data \
  batvi:ubuntu preflight /data
```

## 6) Run Pipeline / 执行主流程

```bash
docker run --rm -it --platform linux/amd64 \
  -v /absolute/path/to/processing:/data \
  batvi:ubuntu run /data -t 8 -l batvi.log
```

## 7) Config Essentials / 配置关键项

In `batviconfig.txt`:

- `BLAST_PATH=/usr/bin`
- `BWA_PATH=/usr/bin`
- `SAMTOOLS_PATH=/usr/bin`
- `BEDTOOLS_PATH=/usr/bin`

Picard is supported via:

- `PICARD_PATH/SamToFastq.jar` (legacy), or
- `picard` CLI, or
- `/usr/share/java/picard.jar`.

Template file: `/Users/kenneth/BatVI/batviconfig.template.txt`

## 8) Troubleshooting / 常见问题

1. `-msse2` compile issues on macOS host:
   - Build/run in Docker `linux/amd64`.
2. Docker runs but preflight fails:
   - Check index paths in `batviconfig.txt`.
   - Ensure `filelist.txt` entries match mounted filenames.
3. Picard not found:
   - In container this is preinstalled.
   - On host, use `picard` CLI or set `PICARD_PATH`.
