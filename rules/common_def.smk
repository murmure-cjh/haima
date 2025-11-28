# =============================================================================
# rules/common.smk - 通用函数和工具模块
# 
# 功能：提供Snakemake流程所需的通用函数，包括：
#   - 命令行参数解析
#   - 目录路径管理
#   - 样本配置管理
#   - 文件路径获取
#   - 流程判断逻辑
# =============================================================================

import os
import glob
import yaml


# =============================================================================
# 全局 Shell 安全设置
# =============================================================================

shell.prefix("set -eo pipefail; ")

# =============================================================================
# 命令行参数处理函数
# =============================================================================

def parse_args():
    """
    解析命令行参数
    
    Returns:
        tuple: (已知参数, 未知参数)
    """
    import argparse
    
    parser = argparse.ArgumentParser(description='Haima分析流程')
    
    # 默认使用最新的yaml文件
    default_config = get_latest_sample_config()
    
    parser.add_argument(
        '--sample-config', 
        type=str, 
        default=default_config,
        help='样本配置文件路径 (默认: 最新的yaml文件)'
    )
    
    return parser.parse_known_args()

# =============================================================================
# 目录路径管理函数
# =============================================================================

def create_sample_dirs(wildcards):
    """
    创建样本所需的所有目录结构
    
    Args:
        wildcards: Snakemake通配符对象，包含sample和sed_id
        
    Returns:
        str: 创建目录的shell命令字符串
    """
    sample_name = wildcards.sample
    sed_id = wildcards.sed_id
    
    # 定义所有需要创建的目录
    directories = [
        f"results/{sed_id}/{sample_name}",
        f"results/{sed_id}/{sample_name}/alignment",
        f"results/{sed_id}/{sample_name}/annotation", 
        f"results/{sed_id}/{sample_name}/gender",
        f"results/{sed_id}/{sample_name}/dipin",
        f"results/{sed_id}/{sample_name}/hg38",
        f"results/{sed_id}/{sample_name}/logs",
        f"results/{sed_id}/{sample_name}/qc",
        f"results/{sed_id}/{sample_name}/qc/tmp",
        f"results/{sed_id}/{sample_name}/sma",
        f"results/{sed_id}/{sample_name}/snp",
        f"results/{sed_id}/{sample_name}/tmp",
        f"results/{sed_id}/{sample_name}/variant"
    ]
    
    # 生成创建目录的shell命令
    mkdir_cmds = [f"mkdir -p {dir}" for dir in directories]
    return " && ".join(mkdir_cmds)


def create_sample_dirs_shell(wildcards):
    """返回创建目录的shell命令字符串"""
    return create_sample_dirs(wildcards)


def get_sample_sed_id(sample_name):
    """
    根据样本名获取对应的sed_id
    
    Args:
        sample_name (str): 样本名称
        
    Returns:
        str: sed_id，如果未找到则返回'unknown_sed'
    """
    return SAMPLES.get(sample_name, {}).get('Sed_ID', 'unknown_sed')


def get_sample_dynamic_path(sample_name, relative_path):
    """
    生成包含正确sed_id的完整路径
    
    Args:
        sample_name (str): 样本名称
        relative_path (str): 相对路径
        
    Returns:
        str: 完整路径
    """
    sed_id = get_sample_sed_id(sample_name)
    return f"results/{sed_id}/{sample_name}/{relative_path}"


def get_sample_base_dir(wildcards):
    """
    获取样本基础目录: results/{sed_id}/{sample}
    
    Args:
        wildcards: Snakemake通配符对象或字典
        
    Returns:
        str: 基础目录路径
    """
    # 处理 wildcards 为字典的情况（在顶层调用时）
    if isinstance(wildcards, dict):
        sample_name = wildcards.get('sample', 'unknown_sample')
        sed_id = wildcards.get('sed_id', 'unknown_sed')
    else:
        sample_name = wildcards.sample
        sed_id = wildcards.sed_id
        
    return f"results/{sed_id}/{sample_name}"


def get_sample_result_dir(wildcards):
    """获取样本结果目录: results/{sed_id}/{sample}/result"""
    return f"{get_sample_base_dir(wildcards)}/result"


def get_sample_log_dir(wildcards):
    """获取样本日志目录: results/{sed_id}/{sample}/logs"""
    return f"{get_sample_base_dir(wildcards)}/logs"


def get_sample_tmp_dir(wildcards):
    """获取样本临时目录: results/{sed_id}/{sample}/tmp"""
    return f"{get_sample_base_dir(wildcards)}/tmp"

def get_sample_fastp_dir(wildcards):
    """获取fastp质控目录: results/{sed_id}/{sample}/fastp"""
    return f"{get_sample_base_dir(wildcards)}/fastp"


def get_sample_alignment_dir(wildcards):
    """获取比对分析目录: results/{sed_id}/{sample}/alignment"""
    return f"{get_sample_base_dir(wildcards)}/alignment"


def get_sample_variant_dir(wildcards):
    """获取变异检测目录: results/{sed_id}/{sample}/variant"""
    return f"{get_sample_base_dir(wildcards)}/variant"


def get_sample_annotation_dir(wildcards):
    """获取注释分析目录: results/{sed_id}/{sample}/annotation"""
    return f"{get_sample_base_dir(wildcards)}/annotation"


def get_sample_qc_dir(wildcards):
    """获取质控分析目录: results/{sed_id}/{sample}/qc"""
    return f"{get_sample_base_dir(wildcards)}/qc"


def get_sample_snp_dir(wildcards):
    """获取SNP分析目录: results/{sed_id}/{sample}/snp"""
    return f"{get_sample_base_dir(wildcards)}/snp"


def get_sample_gender_dir(wildcards):
    """获取CNV分析目录: results/{sed_id}/{sample}/gender"""
    return f"{get_sample_base_dir(wildcards)}/gender"


def get_sample_sma_dir(wildcards):
    """获取SMA分析目录: results/{sed_id}/{sample}/sma"""
    return f"{get_sample_base_dir(wildcards)}/sma"


def get_sample_dipin_dir(wildcards):
    """获取Dipin分析目录: results/{sed_id}/{sample}/dipin"""
    return f"{get_sample_base_dir(wildcards)}/dipin"


def get_sample_hg38_dir(wildcards):
    """获取hg38分析目录: results/{sed_id}/{sample}/hg38"""
    return f"{get_sample_base_dir(wildcards)}/hg38"


# =============================================================================
# 样本配置管理函数
# =============================================================================

def get_latest_sample_config():
    """
    获取sample_info/snakemake_sample_yaml/目录下最新的yaml文件
    
    Returns:
        str: 最新的yaml文件路径
    """
    sample_dir = "sample_info/snakemake_sample_yaml/"
    
    # 确保目录存在
    if not os.path.exists(sample_dir):
        print(f"警告: 样本目录 {sample_dir} 不存在，使用默认配置")
        return 'sample_info/snakemake_sample_yaml/20251031_LH00128_0450_B2335CGLT4_clinical_qc_haima_1.yaml'
    
    # 查找所有的yaml文件
    yaml_files = glob.glob(os.path.join(sample_dir, "*.yaml")) + \
                 glob.glob(os.path.join(sample_dir, "*.yml"))
    
    if not yaml_files:
        print(f"警告: 在 {sample_dir} 中未找到yaml文件，使用默认配置")
        return 'sample_info/snakemake_sample_yaml/20251031_LH00128_0450_B2335CGLT4_clinical_qc_haima_1.yaml'
    
    # 按修改时间排序，获取最新的文件
    latest_file = max(yaml_files, key=os.path.getmtime)
    print(f"使用最新的样本配置文件: {latest_file}")
    
    return latest_file


def load_sample_config(config_path=None):
    """
    加载样本配置
    
    Args:
        config_path (str, optional): 配置文件路径，默认为None使用最新文件
        
    Returns:
        dict: 样本配置字典
        
    Raises:
        Exception: 当配置文件加载失败时抛出
    """
    if config_path is None:
        config_path = get_latest_sample_config()
    
    print(f"正在加载样本配置: {config_path}")
    
    try:
        with open(config_path, 'r') as f:
            samples_config = yaml.safe_load(f)
        
        # 移除顶层的Sed_ID，因为它已经在每个样本中定义了
        if 'Sed_ID' in samples_config and len(samples_config) > 1:
            # 如果Sed_ID是顶层键且还有其他样本，则移除它
            sed_id = samples_config.pop('Sed_ID', None)
            print(f"发现顶层Sed_ID: {sed_id}")
        
        print(f"成功加载 {len(samples_config)} 个样本")
        
        # 验证每个样本的配置
        for sample_id, sample_info in samples_config.items():
            validate_sample_config(sample_id, sample_info)
            log_sample_info(sample_id, sample_info)
            
        return samples_config
        
    except Exception as e:
        print(f"错误: 无法加载样本配置文件 {config_path}: {e}")
        raise


# =============================================================================
# 样本属性获取函数
# =============================================================================


    

def get_sample_attribute(wildcards, attr):
    """
    安全获取样本属性
    
    Args:
        wildcards: Snakemake通配符对象
        attr (str): 属性名称
        
    Returns:
        str: 属性值，如果不存在则返回空字符串
    """
    return SAMPLES.get(wildcards.sample, {}).get(attr, "")


def get_fastq_path(wildcards, read_type):
    """
    获取fastq文件路径（使用相对路径）
    
    Args:
        wildcards: Snakemake通配符对象
        read_type (str): 读取类型，'R1'或'R2'
        
    Returns:
        str or None: fastq文件相对路径，如果不存在则返回None
    """
    sample_info = SAMPLES.get(wildcards.sample, {})
    
    if read_type == 'R1':
        raw_path = sample_info.get('Local_R1') or sample_info.get('raw_path_r1')
    else:
        raw_path = sample_info.get('Local_R2') or sample_info.get('raw_path_r2')
    
    if not raw_path:
        return None
    
    # 确保路径是绝对路径，然后转换为相对路径
    abs_path = os.path.abspath(raw_path)
    return os.path.relpath(abs_path)


def get_panel_type(wildcards):
    """
    获取panel类型
    
    Args:
        wildcards: Snakemake通配符对象
        
    Returns:
        str: panel类型
    """
    return SAMPLES.get(wildcards.sample, {}).get('Panel', '')


def get_sample_id(wildcards):
    """
    获取样本ID
    
    Args:
        wildcards: Snakemake通配符对象
        
    Returns:
        str: 样本ID
    """
    return SAMPLES.get(wildcards.sample, {}).get('Sample_ID', wildcards.sample)


# =============================================================================
# BED文件选择函数
# =============================================================================

def get_bed_file(wildcards):
    """
    动态获取BED文件路径
    
    Args:
        wildcards: Snakemake通配符对象
        
    Returns:
        str: BED文件路径
    """
    panel = get_panel_type(wildcards)
    bed_file = config['BED_FILE_DICT'].get(panel)
    return bed_file


def get_bed_qc_file(wildcards):
    """
    动态获取QC BED文件路径
    
    Args:
        wildcards: Snakemake通配符对象
        
    Returns:
        str: QC BED文件路径
    """
    panel = get_panel_type(wildcards)
    bed_file = config['BED_FILE_QC_DICT'].get(panel)
    return bed_file


# =============================================================================
# 性别分析参数函数
# =============================================================================

def get_gender_params(wildcards):
    """
    根据panel类型获取ngs-bits SampleGender参数
    
    Args:
        wildcards: Snakemake通配符对象
        
    Returns:
        str: ngs-bits SampleGender参数
    """
    panel = get_panel_type(wildcards)
    
    # 从配置文件中获取参数，如果未找到则使用默认值
    gender_params = config.get('GENDER_PARAMS', {})
    
    # 优先查找特定panel的参数，如果没有则使用默认参数
    params = gender_params.get(panel, gender_params.get('default', ''))
    
    return params


# =============================================================================
# 流程判断函数
# =============================================================================

def is_miniWES(wildcards):
    """
    检查样本是否为miniWES流程
    
    Args:
        wildcards: Snakemake通配符对象
        
    Returns:
        bool: 是否为miniWES流程
    """
    panel = get_panel_type(wildcards)
    is_miniwes = panel == 'NAD2403_miniWES'
    print(f"样本 {wildcards.sample} Panel: {panel}, miniWES流程: {is_miniwes}")
    
    return is_miniwes


def get_miniwes_samples():
    """
    获取所有miniWES样本列表
    
    Returns:
        list: miniWES样本名称列表
    """
    return [
        sample for sample, info in SAMPLES.items() 
        if info.get('Panel') == 'NAD2403_miniWES'
    ]


# =============================================================================
# SMA分析相关函数
# =============================================================================

def get_gender_from_file(gender_file):
    """
    从性别文件中获取性别信息（第二行第二列）
    
    Args:
        gender_file (str): 性别文件路径
        
    Returns:
        str: 性别信息 ('male', 'female', 或 'unknown')
    """
    try:
        with open(gender_file, 'r') as f:
            lines = f.readlines()
            if len(lines) >= 2:
                # 第二行第二列
                gender = lines[1].strip().split('\t')[1].lower()
                if gender in ['male', 'female']:
                    return gender
                else:
                    return 'unknown'
            else:
                return 'unknown'
                
    except Exception as e:
        return 'unknown'


def get_dmd_bam_dir(wildcards):
    """
    根据性别获取DMD BAM文件目录
    
    Args:
        wildcards: Snakemake通配符对象
        
    Returns:
        str: DMD BAM文件目录路径
    """
    gender_file = f"{get_sample_result_dir(wildcards)}/CNV/{wildcards.sample}.gender.txt"
    gender = get_gender_from_file(gender_file)
    
    if gender == 'male':
        return config['dmd_bd_man_dir']
    elif gender == 'female':
        return config['dmd_bd_woman_dir']
    else:
        # 如果性别未知，默认使用男性目录
        return config['dmd_bd_man_dir']


def get_dmd_bam_files(wildcards):
    """
    获取DMD目录下的所有BAM文件
    
    Args:
        wildcards: Snakemake通配符对象
        
    Returns:
        list: BAM文件路径列表
        
    Raises:
        ValueError: 当在DMD目录中未找到BAM文件时抛出
    """
    dmd_dir = get_dmd_bam_dir(wildcards)
    import glob
    
    bam_files = glob.glob(f"{dmd_dir}/*.bam")
    
    if not bam_files:
        raise ValueError(f"在DMD目录 {dmd_dir} 中未找到BAM文件")
    
    return bam_files


# =============================================================================
# 工具函数
# =============================================================================

def validate_sample_config(sample_id, sample_info):
    """
    验证样本配置是否完整
    
    Args:
        sample_id (str): 样本ID
        sample_info (dict): 样本信息字典
        
    Returns:
        bool: 验证是否通过
        
    Raises:
        ValueError: 当缺少必要字段时抛出
    """
    required_fields = ['Local_R1', 'Local_R2', 'Panel', 'Sed_ID']
    missing_fields = [field for field in required_fields if field not in sample_info]
    
    if missing_fields:
        raise ValueError(f"样本 {sample_id} 缺少必要字段: {missing_fields}")
    
    # 检查文件是否存在
    for read_type in ['Local_R1', 'Local_R2']:
        file_path = sample_info[read_type]
        if not os.path.exists(file_path):
            print(f"警告: 样本 {sample_id} 的 {read_type} 文件不存在: {file_path}")
    
    return True


def log_sample_info(sample_id, sample_info):
    """
    记录样本信息 - 美化版本
    
    Args:
        sample_id (str): 样本ID
        sample_info (dict): 样本信息字典
    """
    # 颜色代码 (ANSI)
    COLORS = {
        'blue': '\033[94m',
        'green': '\033[92m',
        'yellow': '\033[93m',
        'red': '\033[91m',
        'bold': '\033[1m',
        'end': '\033[0m'
    }
    
    # 获取数据
    panel = sample_info.get('Panel', '未知')
    r1_path = sample_info.get('Local_R1', '未知')
    r2_path = sample_info.get('Local_R2', '未知')
    sed_id = sample_info.get('Sed_ID', '未知')
    is_mini_wes = panel == 'NAD2403_miniWES'
    
    # 构建输出
    header = f"{COLORS['bold']}{COLORS['blue']}📊 样本信息{COLORS['end']}"
    separator = f"{COLORS['blue']}{'='*50}{COLORS['end']}"
    
    print()
    print()
    print(separator)
    print(header)
    print(separator)
    
    print(f"{COLORS['bold']}样本ID:{COLORS['end']} {COLORS['green']}{sample_id}{COLORS['end']}")
    print(f"{COLORS['bold']}Sed_ID:{COLORS['end']} {sed_id}")
    print(f"{COLORS['bold']}Panel类型:{COLORS['end']} {COLORS['yellow']}{panel}{COLORS['end']}")
    
    # 文件路径显示
    print(f"{COLORS['bold']}测序文件:{COLORS['end']}")
    print(f"  ├─ R1: {r1_path}")
    print(f"  └─ R2: {r2_path}")
    
    # miniWES状态
    status_color = COLORS['green'] if is_mini_wes else COLORS['yellow']
    status_text = "✅ 是" if is_mini_wes else "❌ 否"
    print(f"{COLORS['bold']}是否为miniWES:{COLORS['end']} {status_color}{status_text}{COLORS['end']}")
    
    print(separator)
    print()  # 空行分隔