#!/usr/bin/env python3
"""
海马WES数据处理结果修饰脚本
对ACMG分类器输出结果进行格式化和字段调整
支持合并两个JSON文件并根据chr+position去重
输出为TXT表格格式
"""

import json
import sys
import re
from typing import Dict, List, Any

def extract_exon_number(exon_str: str) -> str:
    """提取exon数字"""
    if exon_str in ['.', '-', '']:
        return None
    # 匹配第一个数字
    match = re.search(r'(\d+)', exon_str)
    if match:
        return f"exon{match.group(1)}"
    return None

def extract_intron_number(intron_str: str) -> str:
    """提取intron数字"""
    if intron_str in ['.', '-', '']:
        return None
    # 匹配第一个数字
    match = re.search(r'(\d+)', intron_str)
    if match:
        return f"intron{match.group(1)}"
    return None

def extract_hgvs_content(hgvs_str: str) -> str:
    """提取HGVSc冒号后的内容"""
    if hgvs_str in ['.', '-', '']:
        return hgvs_str
    if ':' in hgvs_str:
        return hgvs_str.split(':', 1)[1]
    return hgvs_str

def extract_aa_change(aa_change_str: str) -> str:
    """提取AAChange最后一个冒号后的内容"""
    if aa_change_str in ['.', '-', '']:
        return '-'
    if ':' in aa_change_str:
        return aa_change_str.split(':')[-1]
    return aa_change_str

def transform_variant(variant: Dict) -> Dict:
    """
    转换单个变异记录的字段
    
    Args:
        variant: 原始变异记录
        
    Returns:
        转换后的变异记录
    """
    transformed = {}
    
    # 按照模板顺序创建所有字段并初始化为空
    template_fields = [
        'disease', 'OMIM', 'Gene', 'Transcript', 'MutName', 'AminoAcidChange', 
        'region', 'zygosity', 'Chr', 'position', 'ref', 'alt', 'rsID', 
        'G1000_AF', 'G1000_EAS_AF', 'dbSNP_AF', 'ESP6500_AF', 'ExAC_AF', 
        'gnomad_AF', 'gnomad_AF_EAS', 'FunctionalChange', 'EnsemblGeneID', 
        'SIFT', 'Polyphen2_HVAR', 'refDepth', 'altDepth', 'Readsratio', 
        'allDepth', 'HGMD2017', 'Clinvar_CLNSIG', 'classification', 
        'HapWES10000_AF', 'Bioinf_classification', 'MutationTaster_pred', 
        'REVEL', 'Spliceai_masked', 'MaxEntScan_Reduced_value'
    ]
    
    # 初始化所有字段为空
    for field in template_fields:
        transformed[field] = ''
    
    # 设置默认值
    transformed['disease'] = 'NN'  # disease列默认值为NN
    transformed['dbSNP_AF'] = '0'  # dbSNP_AF列默认值为0
    transformed['EnsemblGeneID'] = '.'  # EnsemblGeneID列默认值为.
    transformed['Polyphen2_HVAR'] = '.'  # Polyphen2_HVAR列默认值为.
    transformed['REVEL'] = '.' # REVEL列默认值为.
        
    # 填充实际有值的字段
    # 基本信息
    transformed['Chr'] = variant.get('Chr', '')
    transformed['position'] = variant.get('Start', '')
    transformed['ref'] = variant.get('Ref', '')
    transformed['alt'] = variant.get('Alt', '')
    transformed['Gene'] = variant.get('Gene.refGene', '')
    transformed['Transcript'] = variant.get('HGNC_ID', '')
    transformed['MutName'] = extract_hgvs_content(variant.get('HGVSc', ''))
    transformed['AminoAcidChange'] = extract_aa_change(variant.get('AAChange.refGene', ''))
    transformed['rsID'] = variant.get('avsnp151', '')
    
    # 处理region字段 (EXON/INTRON)
    exon_region = extract_exon_number(variant.get('EXON', ''))
    intron_region = extract_intron_number(variant.get('INTRON', ''))
    
    if exon_region:
        transformed['region'] = exon_region
    elif intron_region:
        transformed['region'] = intron_region
    
    # 频率相关字段
    transformed['gnomad_AF'] = variant.get('AF', '')
    transformed['gnomad_AF_EAS'] = variant.get('AF_eas', '')
    transformed['G1000_AF'] = variant.get('1000g2015aug_all', '')
    transformed['G1000_EAS_AF'] = variant.get('1000g2015aug_eas', '')
    transformed['ESP6500_AF'] = variant.get('esp6500siv2_all', '')
    transformed['ExAC_AF'] = variant.get('ExAC_ALL', '')
    
    # 预测相关字段
    transformed['SIFT'] = variant.get('SIFT_pred', '')
    # Polyphen2_HVAR 已经有默认值'.'
    transformed['MutationTaster_pred'] = variant.get('MutationTaster_pred', '')
    transformed['Spliceai_masked'] = variant.get('Spliceai_masked', '')
    transformed['MaxEntScan_Reduced_value'] = variant.get('MaxEntScan_Reduced_value', '')
    
    # 临床分类字段
    transformed['Clinvar_CLNSIG'] = variant.get('CLNSIG', '')
    
    # 深度相关字段
    transformed['refDepth'] = variant.get('refDepth', 0)
    transformed['altDepth'] = variant.get('altDepth', 0)
    transformed['allDepth'] = variant.get('allDepth', 0)
    transformed['Readsratio'] = variant.get('Readsratio', 0)
    
    # 数据库相关字段
    transformed['HGMD2017'] = variant.get('HGMD', '')
    transformed['OMIM'] = variant.get('OMIM', '')
    transformed['HapWES10000_AF'] = variant.get('HapWES10000_AF', '')
    
    # 分类字段
    transformed['classification'] = variant.get('ACMG_Classification', '')
    transformed['FunctionalChange'] = variant.get('mutation_type_cn', '')
    transformed['Bioinf_classification'] = variant.get('classification_cn', '')
    transformed['zygosity'] = variant.get('zygosity_cn', '')
    
    return transformed

def transform_additional_variant(variant):
    """
    处理额外导入的SNP文件记录 (修正版：恢复原始键名读取逻辑)
    """
    # 1. 定义默认值 (如果取不到值，将回退到这些值)
    default_values = {
        'disease': 'NN',
        'OMIM': '.',
        'Gene': '.',
        'Transcript': '.',
        'MutName': '.',
        'AminoAcidChange': '.',
        'region': '.',
        'zygosity': '.',              
        'FunctionalChange': '.',
        'EnsemblGeneID': '.',
        'rsID': '.',
        'SIFT': '.',
        'Polyphen2_HVAR': '.',
        'MutationTaster_pred': '.',
        'REVEL': '.',
        'Spliceai_masked': '.',
        'MaxEntScan_Reduced_value': '.',
        'HGMD2017': '.',
        'Clinvar_CLNSIG': '.',
        'classification': '.',
        'Bioinf_classification': '.',
        'G1000_AF': '0',
        'G1000_EAS_AF': '0',
        'dbSNP_AF': '0',
        'ESP6500_AF': '0',
        'ExAC_AF': '0',
        'gnomad_AF': '0',
        'gnomad_AF_EAS': '0',
        'HapWES10000_AF': '0',
        'refDepth': '0',
        'altDepth': '0',
        'allDepth': '0',
        'Readsratio': '0'
    }

    # 2. 辅助函数：安全获取非空值
    # 逻辑：优先取 key1，没有则取 key2，如果取出的值为空字符串或None，则返回默认值
    def get_val(keys, default_key):
        val = None
        # 如果传入的是单个字符串，转为列表
        if isinstance(keys, str): keys = [keys]
        
        for k in keys:
            if k in variant and variant[k] not in [None, ""]:
                val = variant[k]
                break
        
        # 如果所有key都没取到有效值，返回预设的默认值
        return val if val is not None else default_values.get(default_key, '.')

    transformed = {}

    # --- 核心修复部分：恢复原始脚本的键名映射 ---

    # 杂合性：优先读取 'zygosity_cn' (原始脚本逻辑)，其次尝试 'zygosity'
    transformed['zygosity'] = get_val(['zygosity_cn', 'zygosity'], 'zygosity')

    # 变异功能：优先读取 'mutation_type_cn'
    transformed['FunctionalChange'] = get_val(['mutation_type_cn', 'Func.refGene', 'ExonicFunc.refGene'], 'FunctionalChange')

    # 生信预测分类：优先读取 'classification_cn'
    transformed['Bioinf_classification'] = get_val(['classification_cn', 'Bioinf_classification'], 'Bioinf_classification')

    # ACMG分类：优先读取 'ACMG_Classification'
    transformed['classification'] = get_val(['ACMG_Classification', 'classification'], 'classification')

    # 基因信息：部分文件可能是 Gene.refGene 或 Gene
    transformed['Gene'] = get_val(['Gene.refGene', 'Gene'], 'Gene')
    
    # --- 下面是常规字段的处理 ---

    transformed['sample_no'] = variant.get('sample_no', '') # 既然是处理结果，这里可能需要外部传入或保持原样
    transformed['disease'] = get_val('disease', 'disease')
    transformed['OMIM'] = get_val('OMIM', 'OMIM')
    
    # 转录本处理 (尝试解析)
    transcript = get_val(['GeneDetail.refGene', 'Transcript'], 'Transcript')
    # 如果没拿到，尝试从 AAChange 中解析 (参考原始脚本逻辑)
    if transcript == '.' and 'AAChange.refGene' in variant:
        aa_str = variant['AAChange.refGene']
        if aa_str and ':' in aa_str:
            parts = aa_str.split(':')
            if len(parts) > 1:
                transcript = parts[1]
    transformed['Transcript'] = transcript

    # 变异名称提取
    transformed['MutName'] = extract_hgvs_content(variant.get('HGVSc', '.'))
    transformed['AminoAcidChange'] = extract_aa_change(variant.get('AAChange.refGene', '.'))

    # 坐标
    transformed['Chr'] = variant.get('Chr', '.')
    transformed['position'] = variant.get('Start', '.')
    transformed['ref'] = variant.get('Ref', '.')

    raw_alt = variant.get('Alt', '.')
    if raw_alt == '.':
        transformed['alt'] = transformed['ref']
    else:
        transformed['alt'] = raw_alt
    # 区域判断 (EXON/INTRON) - 参考原始逻辑
    region_val = '.'
    if 'EXON' in variant and variant['EXON']:
        region_val = extract_exon_number(variant['EXON'])
    elif 'INTRON' in variant and variant['INTRON']:
        region_val = extract_intron_number(variant['INTRON'])
    if not region_val or region_val == '.':
         region_val = variant.get('Func.refGene', default_values['region'])
    transformed['region'] = region_val

    # 数据库ID
    transformed['rsID'] = get_val(['avsnp151', 'avsnp147', 'rsID'], 'rsID')
    transformed['HGMD2017'] = get_val(['HGMD', 'HGMD2017'], 'HGMD2017')
    transformed['Clinvar_CLNSIG'] = get_val(['CLNSIG', 'Clinvar_CLNSIG'], 'Clinvar_CLNSIG')
    transformed['EnsemblGeneID'] = get_val('EnsemblGeneID', 'EnsemblGeneID')

    # 频率 (AF)
    transformed['G1000_AF'] = get_val(['1000g2015aug_all', 'G1000_AF'], 'G1000_AF')
    transformed['G1000_EAS_AF'] = get_val(['1000g2015aug_eas', 'G1000_EAS_AF'], 'G1000_EAS_AF')
    transformed['dbSNP_AF'] = get_val('dbSNP_AF', 'dbSNP_AF') # 如果源文件没有这个key，会自动填0
    transformed['ESP6500_AF'] = get_val(['esp6500siv2_all', 'ESP6500_AF'], 'ESP6500_AF')
    transformed['ExAC_AF'] = get_val(['ExAC_ALL', 'ExAC_AF'], 'ExAC_AF')
    transformed['gnomad_AF'] = get_val(['gnomad211_all', 'AF', 'gnomad_AF'], 'gnomad_AF')
    transformed['gnomad_AF_EAS'] = get_val(['gnomad211_eas', 'AF_eas', 'gnomad_AF_EAS'], 'gnomad_AF_EAS')
    transformed['HapWES10000_AF'] = get_val('HapWES10000_AF', 'HapWES10000_AF')

    # 软件预测分数
    transformed['SIFT'] = get_val(['SIFT_pred', 'SIFT_score', 'SIFT'], 'SIFT')
    transformed['Polyphen2_HVAR'] = get_val(['Polyphen2_HVAR_score', 'Polyphen2_HVAR'], 'Polyphen2_HVAR')
    transformed['MutationTaster_pred'] = get_val('MutationTaster_pred', 'MutationTaster_pred')
    transformed['REVEL'] = get_val(['REVEL_score', 'REVEL'], 'REVEL')
    transformed['Spliceai_masked'] = get_val('Spliceai_masked', 'Spliceai_masked')
    transformed['MaxEntScan_Reduced_value'] = get_val('MaxEntScan_Reduced_value', 'MaxEntScan_Reduced_value')

    # 深度 Depth
    transformed['refDepth'] = get_val('refDepth', 'refDepth')
    transformed['altDepth'] = get_val('altDepth', 'altDepth')
    transformed['allDepth'] = get_val('allDepth', 'allDepth')
    transformed['Readsratio'] = get_val('Readsratio', 'Readsratio')

    return transformed

def process_json_files(input_file: str, haima_snp_file: str, output_file: str):
    """
    处理两个JSON文件，转换所有变异记录并合并
    
    Args:
        input_file: 主输入文件路径
        haima_snp_file: 额外输入文件路径
        output_file: 输出文件路径
    """
    print(f"读取主输入文件: {input_file}")
    
    try:
        with open(input_file, 'r', encoding='utf-8') as f:
            main_data = json.load(f)
    except Exception as e:
        print(f"读取主文件错误: {e}")
        return
    
    if 'data' not in main_data:
        print("错误: 主输入文件格式不正确，缺少 'data' 字段")
        return
    
    # 获取样本ID
    sample_id = main_data.get('sample_id', 'Unknown')
    print(f"样本ID: {sample_id}")
    
    # 转换主文件的变异记录
    print("转换主文件变异记录...")
    transformed_data = []
    seen_positions = set()  # 用于跟踪已处理的chr+position组合
    
    for variant in main_data['data']:
        transformed_variant = transform_variant(variant)
        # 创建位置标识符
        pos_key = (transformed_variant['Chr'], transformed_variant['position'])
        if pos_key not in seen_positions:
            transformed_data.append(transformed_variant)
            seen_positions.add(pos_key)
    
    # 如果有额外文件，处理并合并
    if haima_snp_file:
        print(f"读取额外文件: {haima_snp_file}")
        try:
            with open(haima_snp_file, 'r', encoding='utf-8') as f:
                haima_snp_data = json.load(f)
        except Exception as e:
            print(f"读取额外文件错误: {e}")
            return
        
        if 'data' not in haima_snp_data:
            print("错误: 额外文件格式不正确，缺少 'data' 字段")
            return
        
        print("转换额外文件变异记录...")
        haima_snp_count = 0
        for variant in haima_snp_data['data']:
            transformed_variant = transform_additional_variant(variant)
            # 创建位置标识符
            pos_key = (transformed_variant['Chr'], transformed_variant['position'])
            if pos_key not in seen_positions:
                transformed_data.append(transformed_variant)
                seen_positions.add(pos_key)
                haima_snp_count += 1
        
        print(f"从额外文件添加了 {haima_snp_count} 个新变异记录")
    
    # 按照模板顺序定义输出列（添加sample_no作为第一列）
    columns = [
        'sample_no', 'disease', 'OMIM', 'Gene', 'Transcript', 'MutName', 
        'AminoAcidChange', 'region', 'zygosity', 'Chr', 'position', 'ref', 
        'alt', 'rsID', 'G1000_AF', 'G1000_EAS_AF', 'dbSNP_AF', 'ESP6500_AF', 
        'ExAC_AF', 'gnomad_AF', 'gnomad_AF_EAS', 'FunctionalChange', 
        'EnsemblGeneID', 'SIFT', 'Polyphen2_HVAR', 'refDepth', 'altDepth', 
        'Readsratio', 'allDepth', 'HGMD2017', 'Clinvar_CLNSIG', 'classification', 
        'HapWES10000_AF', 'Bioinf_classification', 'MutationTaster_pred', 
        'REVEL', 'Spliceai_masked', 'MaxEntScan_Reduced_value'
    ]
    
    # 写入输出文件（TXT格式）
    print(f"写入输出文件: {output_file}")
    try:
        with open(output_file, 'w', encoding='utf-8') as f:
            # 写入表头
            f.write('\t'.join(columns) + '\n')
            
            # 写入数据行
            for variant in transformed_data:
                row = [sample_id]  # 第一列为样本ID
                for col in columns[1:]:  # 跳过第一列sample_no
                    value = variant.get(col, '')
                    # 将None转换为空字符串，其他值转换为字符串
                    row.append(str(value) if value is not None else '')
                f.write('\t'.join(row) + '\n')
                
        print("处理完成!")
        print(f"总变异数: {len(transformed_data)}")
    except Exception as e:
        print(f"写入文件错误: {e}")

def main():
    """主函数"""
    if len(sys.argv) < 3:
        print("用法: python transform_haima_result.py -i input.json -a haima_snp.json -o output.txt")
        print("示例: python transform_haima_result.py -i haima_output.json -a haima_snp_data.json -o transformed_output.txt")
        sys.exit(1)
    
    # 解析命令行参数
    input_file = None
    haima_snp_file = None
    output_file = None
    
    i = 1
    while i < len(sys.argv):
        if sys.argv[i] == "-i" and i + 1 < len(sys.argv):
            input_file = sys.argv[i + 1]
            i += 2
        elif sys.argv[i] == "-a" and i + 1 < len(sys.argv):
            haima_snp_file = sys.argv[i + 1]
            i += 2
        elif sys.argv[i] == "-o" and i + 1 < len(sys.argv):
            output_file = sys.argv[i + 1]
            i += 2
        else:
            i += 1
    
    if not all([input_file, output_file]):
        print("错误: 必须提供 -i 和 -o 参数")
        sys.exit(1)
    
    # 处理文件
    process_json_files(input_file, haima_snp_file, output_file)

if __name__ == "__main__":
    main()