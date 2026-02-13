# BatVI / BatVI 双语部署文档

## 1) Status / 项目状态

### English

BatVI is maintained here as a **legacy, reproducibility-focused pipeline**.  
For new platform development, use **AccenTrust Veclon**:

- Veclon repo: [VancouverCanada/veclon-platform](https://github.com/VancouverCanada/veclon-platform)
- Positioning: secure orchestration platform across heterogeneous viral integration engines.

### 中文

BatVI 在本仓库定位为**历史流程维护与可复现实验**。  
若要做新平台与长期演进，建议转向 **AccenTrust Veclon**：

- Veclon 仓库：[VancouverCanada/veclon-platform](https://github.com/VancouverCanada/veclon-platform)
- 定位：面向多引擎病毒整合分析的安全编排平台。

## 2) Why Docker First / 为什么优先 Docker

### English

BatVI native modules use legacy x86 compile flags (for example `-msse2`).  
On macOS (especially Apple Silicon), host-native build is unreliable.  
Use the provided Ubuntu container (`linux/amd64`) as the default runtime.

### 中文

BatVI 的原生模块包含旧版 x86 编译参数（如 `-msse2`）。  
在 macOS（尤其 Apple Silicon）上直接本机构建不稳定。  
建议默认使用仓库提供的 Ubuntu 容器（`linux/amd64`）运行。

## 3) Quick Start / 快速开始

### English

1. Build image:

```bash
docker build --platform linux/amd64 -f docker/Dockerfile -t batvi:ubuntu .
```

2. Prepare processing directory:
   - `filelist.txt`
   - `batviconfig.txt` (or use `/Users/kenneth/BatVI/batviconfig.template.txt` as template)
   - input FASTQ files

3. Run preflight then run pipeline:

```bash
docker run --rm -it --platform linux/amd64 \
  -v /absolute/path/to/processing:/data \
  batvi:ubuntu preflight /data

docker run --rm -it --platform linux/amd64 \
  -v /absolute/path/to/processing:/data \
  batvi:ubuntu run /data -t 8 -l batvi.log
```

### 中文

1. 构建镜像：

```bash
docker build --platform linux/amd64 -f docker/Dockerfile -t batvi:ubuntu .
```

2. 准备处理目录：
   - `filelist.txt`
   - `batviconfig.txt`（可基于 `/Users/kenneth/BatVI/batviconfig.template.txt` 修改）
   - 输入 FASTQ 文件

3. 先做预检，再执行流程：

```bash
docker run --rm -it --platform linux/amd64 \
  -v /absolute/path/to/processing:/data \
  batvi:ubuntu preflight /data

docker run --rm -it --platform linux/amd64 \
  -v /absolute/path/to/processing:/data \
  batvi:ubuntu run /data -t 8 -l batvi.log
```

## 4) Minimal Config Notes / 最小配置说明

### English

Inside Docker, use:

- `BLAST_PATH=/usr/bin`
- `BWA_PATH=/usr/bin`
- `SAMTOOLS_PATH=/usr/bin`
- `BEDTOOLS_PATH=/usr/bin`

Picard conversion is compatible with:

- legacy `PICARD_PATH/SamToFastq.jar`, or
- `picard` CLI, or
- `/usr/share/java/picard.jar`.

### 中文

在 Docker 内建议配置：

- `BLAST_PATH=/usr/bin`
- `BWA_PATH=/usr/bin`
- `SAMTOOLS_PATH=/usr/bin`
- `BEDTOOLS_PATH=/usr/bin`

Picard 转换已兼容以下三种模式：

- 旧版 `PICARD_PATH/SamToFastq.jar`
- `picard` 命令
- `/usr/share/java/picard.jar`

## 5) Repository Map / 仓库结构

- `/Users/kenneth/BatVI/docker/Dockerfile`: Ubuntu runtime image
- `/Users/kenneth/BatVI/docker/entrypoint.sh`: container command entrypoint
- `/Users/kenneth/BatVI/docker/README.md`: detailed container deployment guide
- `/Users/kenneth/BatVI/preflight_check.sh`: environment validation
- `/Users/kenneth/BatVI/call_integrations.sh`: main pipeline entry
- `/Users/kenneth/BatVI/picard_samtofastq.sh`: Picard compatibility wrapper

## 6) Legacy References / 旧版参考

- Original manual: `/Users/kenneth/BatVI/README.txt`
- Progress tracker: `/Users/kenneth/BatVI/PROJECT_PROGRESS.md`
