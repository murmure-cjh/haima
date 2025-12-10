# 🧬 Haima NGS Analysis Workflow (Snakemake)
这是一个基于 Snakemake 构建的 NGS 自动化分析流程（针对 Haima 项目）。该流程涵盖了从原始数据下载、预处理、变异检测（Calling）、注释（Annotation）到最终报告生成的全过程。流程特别集成了 SMA、地贫（Dipin）以及药物基因组的分析模块。


## 📋 目录结构
```
.
├── Snakefile                   # Snakemake 主入口文件
├── config.yaml                 # 流程主配置文件（包含参考基因组路径、软件路径等）
├── rules/                      # 模块化的 Snakemake 规则
│   ├── common_def.smk          # 通用函数定义
│   ├── core_analysis.smk       # 核心分析模块 (QC, Alignment, Variant Calling)
│   ├── annotation_analysis.smk # 注释与下游分析模块
│   └── miniwes.smk             # 小全外显子组特定分析 (SMA, Dipin 等)
├── scripts/                    # Python 脚本集
│   ├── analysis_scripts/       # 生物信息分析核心脚本 (注释, ACMG, 结果处理)
│   └── workflow_scripts/       # 流程控制脚本 (数据下载, 监控, 邮件发送)
├── sample_info/                # 样本信息管理
│   ├── raw_haima_csv/          # 原始样本信息表 (CSV)
│   └── snakemake_sample_yaml/  # 生成的流程配置文件 (YAML)
├── docker/                     # Docker 环境构建文件
├── env_yaml/                   # Conda 环境配置文件
└── raw_data/                   # [临时] 原始 FASTQ 数据存放目录
```

## 📜 关键脚本说明

### 📂 workflow_scripts (流程控制)

位于 `scripts/workflow_scripts/`，主要负责流程的调度、监控和辅助功能。

| 脚本名 | 功能描述 |
| :--- | :--- |
| `haima_preprocess.py` | **预处理核心**：将云端的样本信息表（CSV）转换为 Snakemake 所需的 YAML 配置文件，并处理 FASTQ 数据的下载逻辑。 |
| `monitor.py` | **监控主程序**：负责监控新任务，自动触发 Snakemake 分析流程。集成了发送邮件和合并 QC 的功能模块。 |
| `merge_qc.py` | **QC合并**：当批次内所有样本分析完成后，将 QC 结果与样本信息表合并，生成最终的质控报告（作为邮件附件）。 |
| `send_mail.py` | **邮件通知**：分析结束后的邮件发送模块（沿用旧流程逻辑）。 |

### 📂 analysis_scripts (生信分析)

位于 `scripts/analysis_scripts/`，主要涉及生物信息学分析的具体逻辑，部分脚本被打包进 Docker 镜像中。

| 脚本名 | 功能描述 |
| :--- | :--- |
| `anno_caller.py` | **注释整合**：调用 Annovar 和 VEP 进行变异注释，并整合两者的结果。 |
| `haimaresult.py` | **结果过滤与处理**：核心处理脚本，包含复杂的注释结果解析、附加数据库信息添加以及变异过滤逻辑。 |
| `acmg_classifier.py` | **ACMG 评级**：根据 ACMG 规则对变异位点进行自动评级（被 `haimaresult.py` 调用）。 |
| `haima_snp.py` | **药物基因组**：专门处理药物基因组（PGx）相关的位点分析。 |
| `deal_sma_dipin_result.py` | **SMA/地贫格式化**：处理 SMA 和地贫（Dipin）的分析结果，转换为可上传的格式。 |
| `transform_haima_result.py` | **IT 格式转换**：将 `haimaresult.py` 和 `haima_snp.py` 的输出结果转换为业务系统（IT）所需的特定格式。 |
| `qc_fload_analysis.py` | **质控统计**：处理和分析质控（QC）结果。 |
| `deal_sma_workout.py` | **SMA 分析核心**：SMA 分析的具体计算逻辑（沿用旧流程）。 |


## 🛠️ 环境部署

```
### conda环境
conda activate snakemake
source activate /x03_haplox/users/chenjh/miniforge3/envs/snakemake

### Docker 构建

# 1. 更新 Docker 内的脚本 (解决上下文问题)
cp scripts/analysis_scripts/* docker/scripts/

# 2. 构建镜像
cd docker
docker build -t xiaohaima .
```

## 🚀 流程运行 
### 自动化监控运行

- 使用监控脚本自动处理新任务，并发送邮件等。
```
nohup python scripts/workflow_scripts/monitor.py >> logs/monitor.out 2>&1 &
```

### 手动运行任务
#### 步骤 1: 数据预处理 (Preprocessing)

- 从云端获取样本信息表，并转换为 Snakemake 需要的 YAML 配置文件，同时触发数据下载。

```

# 示例：获取样本表
coscli cp cos://sz-hapseq/rawfq/JX_health/.../sample.csv sample_info/raw_haima_csv/

# 运行预处理脚本
python scripts/workflow_scripts/haima_preprocess.py \
    -i sample_info/raw_haima_csv/sample_list.csv \
    -o sample_info/snakemake_sample_yaml/sample_config.yaml

```

### 步骤 2: 流程检查 (Dry Run)

- 在正式运行前，建议检查语法和执行计划。

```
# 语法检查
snakemake -n --lint

# 查看执行计划 (Dry run)
snakemake --configfile config.yaml -n -p
```

### 步骤 3: 正式运行

```
# 手动运行或者Debug
# 如果不指定 --config sample_config，默认读取目录中日期最新的 YAML 文件。
# 使用 36 核心运行
snakemake --cores 36 -p

# 指定特定配置文件运行
snakemake --cores 36 -p --config sample_config="sample_info/snakemake_sample_yaml/target.yaml"

# 后台运行
nohup snakemake --cores 36 -p --config sample_config="sample_info/snakemake_sample_yaml/20251120_LH00348_0494_B235VM2LT4_clinical_qc_haima_1.yaml" --rerun-incomplete >> logs/snakemake.log 2>&1 &
```


## ⚙️ 管理与维护

```
# 日志文件 (logs/)
snakemake.log: 主流程运行日志 (查找报错关键词: AttributeError, Error)

monitor.log: 监控脚本运行日志

monitor.out: 下载与标准输出日志

# 查看和杀死snakemake相关的任务
ps aux | grep snakemake
kill -9 id 
# 运行监控脚本并开始自动分析任务
nohup python /haplox/users/chenjh/haima/snakemake/scripts/workflow_scripts/monitor.py >> /haplox/users/chenjh/haima/snakemake/logs/monitor.out 2>&1 &

# 目录解锁
snakemake --unlock
```
