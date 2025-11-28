#!/usr/bin/env python3
# encoding=utf-8
"""
海马SNP分析脚本

功能概述:
    - 从BAM文件中提取目标区域的SNP信息
    - 使用bcftools进行变异检测和基因型识别
    - 整合注释信息生成标准化的SNP分析结果

主要处理流程:
    1. 读取药物靶点BED文件
    2. 从BAM文件中提取目标区域
    3. 使用bcftools进行变异检测
    4. 解析VCF文件获取基因型信息
    5. 整合注释信息生成最终输出

依赖工具:
    - samtools: BAM文件处理
    - bcftools: 变异检测和VCF文件生成

输入文件:
    - BAM文件: 重比对后的BAM文件
    - BED文件: 药物靶点区域定义
    - 参考基因组: 用于比对和变异检测
    - 注释JSON: 包含变异功能注释信息

输出文件:
    - JSON格式的SNP分析结果
"""

import sys
import os
import json
import subprocess
import re
from collections import defaultdict


# =============================================================================
# 工具函数模块
# =============================================================================

def run_command(cmd, description=""):
    """
    运行系统命令并检查执行状态
    
    参数:
        cmd: 要执行的命令字符串
        description: 命令描述，用于日志输出
        
    返回:
        bool: 成功返回True，失败返回False
    """
    print(f"执行命令: {description}")
    print(f"完整命令: {cmd}")
    
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    
    if result.returncode != 0:
        print(f"命令执行失败 {description}: {result.stderr}")
        return False
    
    return True


def read_bed_file(bed_file_path):
    """
    读取BED文件，返回位置到rsID和gene的映射
    
    BED文件格式:
        chr start end rs_id gene_name
        示例: chr1 1000 1001 rs1234 CYP2D6
    
    参数:
        bed_file_path: BED文件路径
        
    返回:
        dict: 位置信息字典，格式为 {chr:pos: {rs: rs_id, gene: gene_name}}
    """
    drug_target_data = {}
    
    try:
        with open(bed_file_path, 'r') as file:
            for line_num, line in enumerate(file, 1):
                fields = line.strip().split()
                
                # 验证字段数量
                if len(fields) < 4:
                    print(f"警告: 第{line_num}行字段不足，跳过")
                    continue
                
                chrom, start, end, rs_id = fields[0], fields[1], fields[2], fields[3]
                # 获取gene信息（第五列，可选）
                gene_name = fields[4] if len(fields) >= 5 else ""
                
                # 使用染色体和结束位置作为唯一键
                position_key = f"{chrom}:{end}"
                drug_target_data[position_key] = {
                    'rs': rs_id,
                    'gene': gene_name
                }
        
        print(f"从BED文件加载了 {len(drug_target_data)} 个靶点位置")
        return drug_target_data
        
    except Exception as error:
        print(f"读取BED文件 {bed_file_path} 时出错: {error}")
        return {}


# =============================================================================
# BAM文件处理和变异检测模块
# =============================================================================

def extract_target_regions(bam_file_path, bed_file_path, output_dir, sample_name):
    """
    从BAM文件中提取目标区域
    
    处理步骤:
        1. 使用samtools view提取目标区域的SAM文件
        2. 将SAM转换为BAM格式
        3. 对BAM文件进行排序和索引
    
    参数:
        bam_file_path: 输入BAM文件路径
        bed_file_path: 目标区域BED文件路径
        output_dir: 输出目录
        sample_name: 样本名称
        
    返回:
        tuple: (成功状态, 排序后的BAM文件路径)
    """
    # 定义输出文件路径
    sam_output = os.path.join(output_dir, f"{sample_name}.snp.sam")
    bam_output = os.path.join(output_dir, f"{sample_name}.snp.bam")
    sorted_bam_output = os.path.join(output_dir, f"{sample_name}.snp.sorted.bam")
    
    # 步骤1: 提取目标区域为SAM格式
    extract_cmd = f"samtools view -h {bam_file_path} -L {bed_file_path} > {sam_output}"
    if not run_command(extract_cmd, "提取目标区域"):
        return False, None
    
    # 步骤2: 转换SAM为BAM格式
    convert_cmd = f"samtools view -S -b -o {bam_output} {sam_output}"
    if not run_command(convert_cmd, "SAM转BAM格式"):
        return False, None
    
    # 步骤3: 排序BAM文件
    sort_cmd = f"samtools sort -o {sorted_bam_output} {bam_output}"
    if not run_command(sort_cmd, "排序BAM文件"):
        return False, None
    
    # 步骤4: 创建BAM索引
    index_cmd = f"samtools index {sorted_bam_output}"
    if not run_command(index_cmd, "创建BAM索引"):
        return False, None
    
    return True, sorted_bam_output


def run_bcftools_variant_calling(sorted_bam_path, bed_file_path, reference_genome, output_dir, sample_name):
    """
    使用bcftools进行变异检测
    
    处理步骤:
        1. bcftools mpileup生成pileup文件
        2. bcftools call进行变异检测
        3. 解压生成的VCF文件
    
    参数:
        sorted_bam_path: 排序后的BAM文件路径
        bed_file_path: 目标区域BED文件路径
        reference_genome: 参考基因组文件路径
        output_dir: 输出目录
        sample_name: 样本名称
        
    返回:
        tuple: (成功状态, VCF文件路径)
    """
    vcf_output = os.path.join(output_dir, f"{sample_name}.snp.vcf")
    vcf_gz_output = f"{vcf_output}.gz"
    
    # 构建bcftools命令
    mpileup_cmd = f"bcftools mpileup -f {reference_genome} -R {bed_file_path} {sorted_bam_path}"
    call_cmd = f"bcftools call -m -Oz -o {vcf_gz_output}"
    
    try:
        print("运行bcftools进行变异检测...")
        
        # 使用管道连接mpileup和call命令
        mpileup_process = subprocess.Popen(mpileup_cmd.split(), stdout=subprocess.PIPE)
        call_process = subprocess.Popen(call_cmd.split(), stdin=mpileup_process.stdout)
        
        # 关闭mpileup的输出流，避免阻塞
        mpileup_process.stdout.close()
        
        # 等待call命令完成
        call_process.communicate()
        
        if call_process.returncode != 0:
            print("bcftools变异检测失败")
            return False, None
            
    except Exception as error:
        print(f"运行bcftools时出错: {error}")
        return False, None
    
    # 解压VCF文件
    if not run_command(f"gunzip -f {vcf_gz_output}", "解压VCF文件"):
        return False, None
    
    return True, vcf_output


def process_with_bcftools(sample_name, bam_file_path, bed_file_path, reference_genome, output_dir):
    """
    使用bcftools处理BAM文件，生成VCF并进行基因型识别
    
    参数:
        sample_name: 样本名称
        bam_file_path: BAM文件路径
        bed_file_path: 目标区域BED文件路径
        reference_genome: 参考基因组文件路径
        output_dir: 输出目录
        
    返回:
        str: VCF文件路径，失败返回None
    """
    # 提取目标区域
    success, sorted_bam_path = extract_target_regions(bam_file_path, bed_file_path, output_dir, sample_name)
    if not success:
        return None
    
    # 运行变异检测
    success, vcf_file_path = run_bcftools_variant_calling(
        sorted_bam_path, bed_file_path, reference_genome, output_dir, sample_name
    )
    
    return vcf_file_path if success else None


# =============================================================================
# VCF文件解析模块
# =============================================================================

def parse_genotype_from_vcf_fields(format_fields, sample_fields):
    """
    从VCF格式字段和样本字段中解析基因型信息
    
    参数:
        format_fields: FORMAT字段列表 (如: ['GT', 'AD', 'DP'])
        sample_fields: 样本数据字段列表 (如: ['0/1', '10,5', '15'])
        
    返回:
        dict: 基因型相关信息
    """
    genotype_info = {
        'genotype': "./.",
        'allele_depth': "0,0",
        'depth': "0"
    }
    
    # 查找各字段在FORMAT中的位置
    gt_index = format_fields.index('GT') if 'GT' in format_fields else -1
    ad_index = format_fields.index('AD') if 'AD' in format_fields else -1
    dp_index = format_fields.index('DP') if 'DP' in format_fields else -1
    
    # 提取基因型信息
    if gt_index >= 0 and gt_index < len(sample_fields):
        genotype_info['genotype'] = sample_fields[gt_index]
    
    if ad_index >= 0 and ad_index < len(sample_fields):
        genotype_info['allele_depth'] = sample_fields[ad_index]
    
    if dp_index >= 0 and dp_index < len(sample_fields):
        genotype_info['depth'] = sample_fields[dp_index]
    
    return genotype_info


def determine_zygosity_and_counts(genotype, allele_depth, depth):
    """
    根据基因型确定纯合/杂合状态和等位基因计数
    
    参数:
        genotype: 基因型字符串 (如: '0/1', '1/1')
        allele_depth: 等位基因深度字符串 (如: '10,5')
        depth: 总深度字符串
        
    返回:
        tuple: (纯合状态, 参考等位基因计数, 变异等位基因计数)
    """
    # 纯合参考基因型
    if genotype in ['0/0', '0|0']:
        zygosity = "纯合"
        ref_count = int(depth) if depth.isdigit() else 0
        alt_count = 0
    
    # 纯合变异基因型
    elif genotype in ['1/1', '1|1']:
        zygosity = "纯合"
        ref_count = 0
        alt_count = int(depth) if depth.isdigit() else 0
    
    # 杂合基因型
    elif genotype in ['0/1', '0|1', '1/0', '1|0']:
        zygosity = "杂合"
        ad_values = allele_depth.split(',')
        
        if len(ad_values) >= 2:
            ref_count = int(ad_values[0]) if ad_values[0].isdigit() else 0
            alt_count = int(ad_values[1]) if ad_values[1].isdigit() else 0
        else:
            ref_count = 0
            alt_count = 0
    
    # 其他未知基因型，默认设为纯合
    else:
        zygosity = "纯合"
        ref_count = 0
        alt_count = 0
    
    return zygosity, ref_count, alt_count


def parse_vcf_genotypes(vcf_file_path, drug_target_data):
    """
    解析VCF文件，获取基因型信息
    
    参数:
        vcf_file_path: VCF文件路径
        drug_target_data: 药物靶点数据字典
        
    返回:
        dict: 基因型数据，格式为 {位置: 基因型信息}
    """
    genotype_data = {}
    processed_variants = 0
    
    try:
        with open(vcf_file_path, 'r') as file:
            for line in file:
                # 跳过注释行
                if line.startswith('#'):
                    continue
                
                fields = line.strip().split('\t')
                # 验证VCF行格式
                if len(fields) < 10:
                    continue
                
                # 解析VCF字段
                chrom, pos, variant_id, ref_allele, alt_allele = fields[0:5]
                qual, filter_status, info, format_field = fields[5:9]
                sample_data = fields[9]
                
                # 构建位置键
                position_key = f"{chrom}:{pos}"
                
                # 检查是否在药物靶点中
                if position_key not in drug_target_data:
                    continue
                
                # 解析格式字段和样本数据
                format_fields = format_field.split(':')
                sample_fields = sample_data.split(':')
                
                # 获取基因型信息
                gt_info = parse_genotype_from_vcf_fields(format_fields, sample_fields)
                
                # 确定纯合状态和等位基因计数
                zygosity, ref_count, alt_count = determine_zygosity_and_counts(
                    gt_info['genotype'], gt_info['allele_depth'], gt_info['depth']
                )
                
                # 存储基因型信息
                genotype_data[position_key] = {
                    'chr': chrom,
                    'pos': pos,
                    'ref': ref_allele,
                    'alt': alt_allele,
                    'rsid': drug_target_data.get(position_key, {}).get('rs', ''),
                    'gene': drug_target_data.get(position_key, {}).get('gene', ''),
                    'genotype': gt_info['genotype'],
                    'zygosity': zygosity,
                    'depth': gt_info['depth'],
                    'ref_count': ref_count,
                    'alt_count': alt_count,
                    'vcf_id': variant_id
                }
                
                processed_variants += 1
        
        print(f"从VCF文件中解析了 {processed_variants} 个变异位点")
        return genotype_data
    
    except Exception as error:
        print(f"解析VCF文件 {vcf_file_path} 时出错: {error}")
        return {}


# =============================================================================
# 注释数据处理模块
# =============================================================================

def load_annotation_json(annotation_file_path):
    """
    加载注释JSON文件
    
    注释JSON格式:
        {
            "data": [
                {
                    "Chr": "chr1",
                    "Start": "1000",
                    "HGVSc": "c.123A>T",
                    "AAChange.refGene": "p.Lys42Asn"
                },
                ...
            ]
        }
    
    参数:
        annotation_file_path: 注释JSON文件路径
        
    返回:
        dict: 注释数据字典，格式为 {chr:pos: {HGVSc: value, AAChange: value}}
    """
    annotation_data = {}
    
    try:
        with open(annotation_file_path, 'r') as file:
            data = json.load(file)
        
        # 处理注释数据
        if 'data' in data and isinstance(data['data'], list):
            for item in data['data']:
                chrom = item.get('Chr', '')
                start_pos = item.get('Start', '')
                
                if chrom and start_pos:
                    position_key = f"{chrom}:{start_pos}"
                    annotation_data[position_key] = {
                        'HGVSc': item.get('HGVSc', 'c.0>0'),
                        'AAChange': item.get('AAChange.refGene', 'p.0>0')
                    }
        
        print(f"从JSON文件加载了 {len(annotation_data)} 条注释信息")
        return annotation_data
    
    except Exception as error:
        print(f"加载注释JSON文件 {annotation_file_path} 时出错: {error}")
        return {}


# =============================================================================
# 结果生成模块
# =============================================================================

def calculate_reads_ratio(ref_count, alt_count, total_depth):
    """
    计算变异等位基因读取比例
    
    参数:
        ref_count: 参考等位基因计数
        alt_count: 变异等位基因计数
        total_depth: 总深度
        
    返回:
        float: 变异等位基因比例 (0-1之间)
    """
    try:
        total_depth_int = int(total_depth)
        if total_depth_int > 0:
            return alt_count / total_depth_int
        else:
            return 0.0
    except (ValueError, TypeError):
        return 0.0


def create_default_position_data(chrom, pos, drug_target_info):
    """
    为VCF中未检测到的位点创建默认数据
    
    参数:
        chrom: 染色体
        pos: 位置
        drug_target_info: 药物靶点信息
        
    返回:
        dict: 默认位置数据
    """
    return {
        "Chr": chrom,
        "Start": int(pos),
        "Ref": "N",
        "Alt": "N",
        "Gene": drug_target_info.get('gene', ''),
        "HGVSc": "c.0>0",
        "AAChange.refGene": "p.0>0",
        "avsnp151": drug_target_info.get('rs', ''),
        "GT": "./.",
        "DP": "0",
        "refDepth": 0,
        "altDepth": 0,
        "allDepth": "0",
        "Readsratio": 0,
        "zygosity_cn": "纯合"
    }


def create_detected_position_data(genotype_info, annotation_info, drug_target_info):
    """
    为VCF中检测到的位点创建数据
    
    参数:
        genotype_info: 基因型信息
        annotation_info: 注释信息
        drug_target_info: 药物靶点信息
        
    返回:
        dict: 检测到的位置数据
    """
    ref_count = genotype_info.get('ref_count', 0)
    alt_count = genotype_info.get('alt_count', 0)
    total_depth = genotype_info.get('depth', '0')
    
    # 计算读取比例
    reads_ratio = calculate_reads_ratio(ref_count, alt_count, total_depth)
    
    return {
        "Chr": genotype_info.get('chr', ''),
        "Start": int(genotype_info.get('pos', 0)),
        "Ref": genotype_info.get('ref', 'N'),
        "Alt": genotype_info.get('alt', 'N'),
        "Gene": genotype_info.get('gene', ''),
        "HGVSc": annotation_info.get('HGVSc', 'c.0>0'),
        "AAChange.refGene": annotation_info.get('AAChange', 'p.0>0'),
        "avsnp151": drug_target_info.get('rs', ''),
        "GT": genotype_info.get('genotype', './.'),
        "DP": total_depth,
        "refDepth": ref_count,
        "altDepth": alt_count,
        "allDepth": total_depth,
        "Readsratio": reads_ratio,
        "zygosity_cn": genotype_info.get('zygosity', '纯合')
    }


def generate_output(sample_name, drug_target_data, genotype_data, annotation_data, output_dir):
    """
    生成输出JSON文件，按照指定格式
    
    参数:
        sample_name: 样本名称
        drug_target_data: 药物靶点数据
        genotype_data: 基因型数据
        annotation_data: 注释数据
        output_dir: 输出目录
        
    返回:
        bool: 成功返回True，失败返回False
    """
    output_data = {'data': []}
    
    # 遍历BED文件中的所有位点，确保每个位点都有输出
    for position_key in drug_target_data.keys():
        # 解析染色体和位置
        chrom_pos = position_key.split(':')
        if len(chrom_pos) != 2:
            continue
            
        chrom, pos = chrom_pos
        drug_target_info = drug_target_data.get(position_key, {})
        
        # 获取基因型和注释信息
        genotype_info = genotype_data.get(position_key, {})
        annotation_info = annotation_data.get(position_key, {})
        
        # 根据是否检测到变异创建相应的数据
        if not genotype_info:
            # 未检测到变异，使用默认数据
            position_data = create_default_position_data(chrom, pos, drug_target_info)
        else:
            # 检测到变异，使用实际数据
            position_data = create_detected_position_data(
                genotype_info, annotation_info, drug_target_info
            )
        
        output_data['data'].append(position_data)
    
    # 写入输出JSON文件
    output_file_path = os.path.join(output_dir, f"{sample_name}_snp_results.json")
    
    try:
        with open(output_file_path, 'w', encoding='utf-8') as file:
            json.dump(output_data, file, indent=2, ensure_ascii=False)
        
        print(f"结果已写入: {output_file_path}")
        print(f"输出文件中包含 {len(output_data['data'])} 个位点")
        return True
        
    except Exception as error:
        print(f"写入输出文件 {output_file_path} 时出错: {error}")
        return False


# =============================================================================
# 主函数
# =============================================================================

def main():
    """
    主函数 - 协调整个SNP分析流程
    
    处理流程:
        1. 参数验证和目录创建
        2. 读取药物靶点BED文件
        3. 使用bcftools处理BAM文件
        4. 解析VCF文件获取基因型
        5. 加载注释信息
        6. 生成最终输出
    """
    # 参数验证
    if len(sys.argv) != 7:
        print("""
海马SNP分析工具

用法: 
    python haima_snp.py <样本名称> <BAM文件> <药物BED文件> <参考基因组> <输出目录> <注释JSON文件>

参数说明:
    <样本名称>       样本标识符
    <BAM文件>        重比对后的BAM文件路径
    <药物BED文件>    药物靶点区域定义文件 (BED格式)
    <参考基因组>     参考基因组文件路径 (如: hg19.fa)
    <输出目录>       结果输出目录
    <注释JSON文件>   变异功能注释文件 (JSON格式)

使用示例:
    python haima_snp.py sample001 /path/to/sample.realigned.bam /path/to/drug_targets.bed /path/to/hg19.fa /path/to/output /path/to/annotations.json
        """)
        sys.exit(1)

    # 解析命令行参数
    sample_name, bam_file_path, bed_file_path, reference_genome, output_dir, annotation_file = sys.argv[1:7]
    
    # 创建输出目录
    os.makedirs(output_dir, exist_ok=True)
    
    print("=" * 60)
    print("海马SNP分析流程启动")
    print("=" * 60)
    print(f"样本名称: {sample_name}")
    print(f"BAM文件: {bam_file_path}")
    print(f"药物BED文件: {bed_file_path}")
    print(f"参考基因组: {reference_genome}")
    print(f"输出目录: {output_dir}")
    print(f"注释文件: {annotation_file}")
    print("-" * 60)
    
    # 步骤1: 读取BED文件
    print("步骤1: 读取药物靶点BED文件...")
    drug_target_data = read_bed_file(bed_file_path)
    if not drug_target_data:
        print("错误: 无法从BED文件加载数据")
        sys.exit(1)
    
    # 步骤2: 使用bcftools处理BAM文件
    print("步骤2: 使用bcftools进行变异检测...")
    vcf_file_path = process_with_bcftools(sample_name, bam_file_path, bed_file_path, reference_genome, output_dir)
    if not vcf_file_path:
        print("错误: bcftools处理失败")
        sys.exit(1)
    
    # 步骤3: 解析VCF文件获取基因型
    print("步骤3: 解析VCF文件获取基因型信息...")
    genotype_data = parse_vcf_genotypes(vcf_file_path, drug_target_data)
    
    # 步骤4: 加载注释JSON
    print("步骤4: 加载注释信息...")
    annotation_data = load_annotation_json(annotation_file)
    
    # 步骤5: 生成输出JSON
    print("步骤5: 生成最终输出...")
    if generate_output(sample_name, drug_target_data, genotype_data, annotation_data, output_dir):
        print("=" * 60)
        print("SNP分析完成!")
        print("=" * 60)
    else:
        print("错误: 生成输出失败")
        sys.exit(1)


if __name__ == "__main__":
    main()