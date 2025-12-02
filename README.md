# snakemake.log文件快速查找报错
AttributeError

# 查看和杀死snakemake相关的任务
ps aux | grep snakemake
kill -9 id 
# 运行监控脚本并开始自动分析任务
nohup python /haplox/users/chenjh/haima/snakemake/scripts/workflow_scripts/monitor.py >> /haplox/users/chenjh/haima/snakemake/logs/monitor.out 2>&1 &


# 下面相当于monitor.py中的具体运行命令，也是本地运行snakemake流程的命令
# 激活环境
conda activate snakemake
source activate /haplox/users/chenjh/miniforge3/envs/snakemake
# 获取样本信息表
coscli cp cos://sz-hapseq/rawfq/JX_health/Nova/sample_info_merge/20251120_LH00348_0494_B235VM2LT4_clinical_qc_haima_1.csv  sample_info/raw_haima_csv/

# 数据预处理脚本，包括下载数据(默认下载到raw_data文件夹下)等，从样本信息表中下载fq文件+获取流程配置文件
python scripts/workflow_scripts/haima_preprocess.py -i sample_info/raw_haima_csv/20251120_LH00348_0494_B235VM2LT4_clinical_qc_haima_1.csv -o sample_info/snakemake_sample_yaml/20251120_LH00348_0494_B235VM2LT4_clinical_qc_haima_1.yaml

# 运行流程，使用对应的样本配置文件，如果不指定，默认使用日期最新的配置文件（03总共72个cpu）
snakemake --cores 36 -p
nohup snakemake --cores 36 -p --rerun-incomplete >> logs/snakemake.log 2>&1 &
nohup snakemake --cores 36 -p --config sample_config="sample_info/snakemake_sample_yaml/20251120_LH00348_0494_B235VM2LT4_clinical_qc_haima_1.yaml" --rerun-incomplete >> logs/snakemake.log 2>&1 &


# 当前路径文件夹说明
-rw-rw-r-- 1 chenjh chenjh    4424 11月 21 17:01 config.yaml  整个流程的配置文件，包含注释前步骤依赖的固定路径，如bed文件，参考基因组等，文件中有具体的注释说明
-rw-r--r-- 1 chenjh chenjh   35386 11月 24 08:46 coscli.log   coscli 下载的日志，目前coscli版本不支持自定义日志路径
drwxrwxr-x 3 chenjh chenjh      49 11月 17 17:38 docker       Dockerfile 以及需要导入Docker的脚本
drwxrwxr-x 2 chenjh chenjh      39 11月 18 11:26 env_yaml     conda环境配置文件，目前只有sentieon步骤实际需要这个conda环境，但是每个rule，snakemake都建议添加conda配置，故所有rule都使用这个xiaohaima_x03.yaml这个配置文件（来源于旧流程主要运行环境/x03_haplox/users/wangx/anaconda3/envs/xiaohaima/）
drwxrwxr-x 2 chenjh chenjh   12288 11月 24 08:46 raw_data     所有的下载fq文件，需要定期清理
-rw-rw-r-- 1 chenjh chenjh    2028 11月 21 18:12 README.md    说明文件
drwxrwxr-x 6 chenjh chenjh    4096 11月 24 08:56 results      结果文件路径
drwxrwxr-x 2 chenjh chenjh     119 11月  7 15:36 rules        子snakemake模块 文件有注释 以下为简单说明 common_def.smk:流程运行需要的函数，可能有冗余 core_analysis.smk：核心rule call变异部分的rule 为每个流程都需要运行模块 annotation_analysis.smk：注释和分析流程模块，对calling的注释和信息处理，也为每个流程都需要运行模块 miniwes.smk：小全外需要单独运行的模块，包含hg38基因组版本的call变异，SMA和Dipin分析等模块
drwxrwxr-x 4 chenjh chenjh      66 11月 19 10:37 sample_info  样本信息文件夹，包含raw_haima_csv 从云上下载的样本信息表 snakemake_sample_yaml经过haima_preprocess.py脚本转换后的yaml文件
drwxrwxr-x 2 chenjh chenjh     78 11月 26 18:22 logs          包含了snakemake.log 流程日志，monitor.log监控日志，monitor.out下载日志
drwxrwxr-x 3 chenjh chenjh    4096 11月 21 14:35 scripts      所有脚本的路径，且全部打包进了Docker，更新脚本需要负责脚本到docker/scripts/下，然后更新Docker
-rw-rw-r-- 1 chenjh chenjh    1588 11月 21 11:34 Snakefile    snakemake主流程文件



# 脚本文件夹说明
# /haplox/users/chenjh/haima/snakemake/scripts/analysis_scripts
-rw-rw-r-- 1 chenjh chenjh 20396 11月  6 11:09 acmg_classifier.py        acmg评级规则脚本，导入haimaresult.py脚本中使用
-rw-rw-r-- 1 chenjh chenjh 31743 11月 21 18:01 anno_caller.py            注释脚本，包含annovar与vep注释，并整合结果
-rw-rw-r-- 1 chenjh chenjh 13332 11月 21 14:35 deal_sma_dipin_result.py  处理地贫和sma结果为可上传的格式
-rw-rw-r-- 1 chenjh chenjh  4828 11月 13 09:49 deal_sma_workout.py       SMA分析需要的脚本，来源于旧流程，未改动
-rw-rw-r-- 1 chenjh chenjh 19519 11月 24 10:42 haima_preprocess.py       流程前处理脚本，将云上的样本信息表，处理流程需要的yaml文件和下载fq文件
-rw-rw-r-- 1 chenjh chenjh 23325 11月 19 17:30 haimaresult.py            注释结果处理脚本，较复杂，包含附加数据添加，过滤等步骤
-rw-rw-r-- 1 chenjh chenjh 22907 11月 21 18:08 haima_snp.py              药物基因组脚本
drwxrwxr-x 2 chenjh chenjh   135 11月 11 11:08 __pycache__
-rw-rw-r-- 1 chenjh chenjh 24354 11月 21 18:12 qc_fload_analysis.py      质控结果处理脚本
-rw-rw-r-- 1 chenjh chenjh 13200 11月 21 17:22 transform_haima_result.py 海码结果处理脚本，处理haimaresult.py和haima_snp.py的结果为IT需要的格式


# /haplox/users/chenjh/haima/snakemake/scripts/workflow_scripts

-rwxrwxr-x 1 chenjh chenjh 19519 11月 24 10:42 haima_preprocess.py        预处理脚本，是将拆分给的csv表转换为snakemake流程需要yaml文件
-rw-rw-r-- 1 chenjh chenjh  3899 11月 28 16:20 merge_qc.py                合并质控结果的脚本，在全部样本分析完成后会将所有样本的qc结果和样本csv表合并，然后作为结束邮件的附件
-rwxrwxr-x 1 chenjh chenjh 11373 11月 28 16:20 monitor.py                 监控脚本，send_mail和merge_qc都为里面的模块，负责自动开始分析snakemake流程
drwxrwxr-x 2 chenjh chenjh    81 11月 28 16:21 __pycache__
-rwxrwxr-x 1 chenjh chenjh  5579 11月 28 11:27 send_mail.py               发送邮件脚本，完全来自于旧流程




# 搭建流程经常用到的命令

# 创建更新docker
cp scripts/analysis_scripts/* docker/scripts/  # 因为构建Docker时上下文的问题，所以将要转移的脚本复杂过去
cd docker
docker build -t xiaohaima .

# 激活环境
conda activate snakemake

# 步骤1: 语法检查
snakemake -n --lint

# 步骤2: 检查配置
snakemake --configfile config.yaml -n
snakemake --configfile config.yaml --rerun-incomplete -n 忽略中间文件
# 步骤3: 干运行检查命令
snakemake -n -p
snakemake -n --reason --printshellcmds

# 具体运行
snakemake --cores 36 -p  --config sample_config="指定路径/配置文件.yaml" # 手动指定预处理脚本的结果文件，默认参数为sample_info/snakemake_sample_y
aml/目录下日期最新的yaml文件
snakemake --cores 36 -p --rerun-incomplete    # 从头开始跑
nohup snakemake --cores 36 -p --rerun-incomplete >> logs/snakemake.log 2>&1 &

# 目录解锁
snakemake --unlock

