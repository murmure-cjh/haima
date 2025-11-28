import os
from pathlib import Path
import yaml

# 包含通用函数（包含样本配置管理函数）
include: "rules/common_def.smk"

# 获取命令行参数
args, unknown = parse_args()

# 从config获取样本配置文件路径，如果没有指定则使用默认
sample_config_path = config.get("sample_config", get_latest_sample_config())

# 加载样本配置
SAMPLES = load_sample_config(sample_config_path)

# 包含主配置
configfile: "config.yaml"

# 包含核心分析模块
include: "rules/core_analysis.smk"

# 包含注释分析模块
include: "rules/annotation_analysis.smk"


# 目标规则：生成最终的filter.vcf文件
rule all:
    input:
        [f"results/{SAMPLES[sample]['Sed_ID']}/{sample}/variant/{sample}.filter.vcf" for sample in SAMPLES.keys()],
        [f"results/{SAMPLES[sample]['Sed_ID']}/{sample}/qc/{sample}.sample_information.txt" for sample in SAMPLES.keys()],
        [f"results/{SAMPLES[sample]['Sed_ID']}/{sample}/gender/{sample}.gender.txt" for sample in SAMPLES.keys()],
        [f"results/{SAMPLES[sample]['Sed_ID']}/{sample}/annotation/{sample}_anno.json" for sample in SAMPLES.keys()],
        [f"results/{SAMPLES[sample]['Sed_ID']}/{sample}/annotation/{sample}_haima.json" for sample in SAMPLES.keys()],
        [f"results/{SAMPLES[sample]['Sed_ID']}/{sample}/snp/{sample}_snp_results.json" for sample in SAMPLES.keys()],
        [f"results/{SAMPLES[sample]['Sed_ID']}/{sample}/annotation/{sample}_germline.vcf" for sample in SAMPLES.keys()],
        [f"results/{SAMPLES[sample]['Sed_ID']}/{sample}/upload_success.txt" for sample in SAMPLES.keys()],

        # SMA和DIPIN结果文件的依赖
        [f"results/{SAMPLES[sample]['Sed_ID']}/{sample}/sma/{sample}_SMA_result.csv" for sample in get_miniwes_samples()],
        [f"results/{SAMPLES[sample]['Sed_ID']}/{sample}/dipin/{sample}_dipin_result.csv" for sample in get_miniwes_samples()],
        [f"results/{SAMPLES[sample]['Sed_ID']}/{sample}/upload_miniwes_success.txt" for sample in get_miniwes_samples()],

# 包含miniWES分析模块（如果有miniWES样本）
if get_miniwes_samples():
    include: "rules/miniwes.smk"
    