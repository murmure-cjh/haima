docker build -t miniforge .



snakemake/
├── Snakefile              # 主工作流文件
├── config.yaml           # 配置文件
├── samples.yaml          # 样本配置文件
├── env_yaml/             # 环境文件目录
│   └── xiaohaima_x03.yaml
└── rules/                # 规则和函数模块
    └── common.smk        # 辅助函数模块


# 验证步骤



conda activate snakemake

# 步骤1: 语法检查
snakemake -n --lint

# 步骤2: 检查配置
snakemake --configfile config.yaml -n
snakemake --configfile config.yaml --rerun-incomplete -n 忽略中间文件
# 步骤3: 干运行检查命令
snakemake -n -p
snakemake -n --reason --printshellcmds
# 步骤4: 生成依赖图
snakemake --dag | dot -Tpng > dag.png

# 3. 详细运行（显示执行的命令）
snakemake --cores 16 -p --rerun-incomplete --report report.html# 重新运行生成这些不完整文件的规则

snakemake --cores 36 -p --rerun-incomplete 
nohup snakemake --cores 36 -p --rerun-incomplete > snakemake.log 2>&1 &
snakemake --cores 8 -j 3
# 目录解锁
snakemake --unlock


# 或者强制重新创建环境
snakemake --use-conda --conda-cleanup-pkgs --conda-create-envs-only -n
# 强制重新创建环境
snakemake --use-conda --conda-cleanup-envs


# snakemake conda 环境
snakemake --use-conda --list-conda-envs
# 使用 --conda-create-envs-only 只创建环境而不运行工作流
snakemake --use-conda --conda-create-envs-only --list-conda-envs