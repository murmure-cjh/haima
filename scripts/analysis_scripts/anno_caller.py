#!/usr/bin/env python3
"""
整合的ANNOVAR和VEP注释流程脚本

功能概述:
    - 对输入的VCF文件进行ANNOVAR注释
    - 过滤不需要的染色体变异
    - 进行VEP功能注释
    - 整合两种注释结果并输出TXT和JSON格式

主要处理流程:
    1. ANNOVAR注释 → 染色体过滤 → VEP注释 → 结果整合
    2. 支持多种基因组数据库和插件
    3. 提供灵活的列筛选和格式输出选项

依赖环境:
    - ANNOVAR: 用于基础变异注释
    - VEP: 用于功能效应预测
    - Python包: pandas, numpy, argparse

作者: chenjh
日期: 1.0
"""

import os
import sys
import subprocess
import pandas as pd
import argparse
import re
import numpy as np
import datetime
import json
import yaml
from pathlib import Path


# =============================================================================
# 配置管理 - 从YAML文件加载配置
# =============================================================================

def load_config(config_path):
    """
    从YAML配置文件加载配置
    
    参数:
        config_path: YAML配置文件路径
        
    返回:
        dict: 配置字典
    """
    try:
        with open(config_path, 'r', encoding='utf-8') as f:
            config = yaml.safe_load(f)
        print(f"成功加载配置文件: {config_path}")
        return config
    except Exception as e:
        print(f"错误: 无法加载配置文件 {config_path}: {e}")
        sys.exit(1)


def get_config_value(config, key_path, default=None):
    """
    安全获取嵌套配置值
    
    参数:
        config: 配置字典
        key_path: 配置键路径，用点分隔，如 'annotation.annovar_script'
        default: 默认值
        
    返回:
        any: 配置值或默认值
    """
    keys = key_path.split('.')
    current = config
    
    for key in keys:
        if isinstance(current, dict) and key in current:
            current = current[key]
        else:
            print(f"警告: 配置项 {key_path} 不存在，使用默认值: {default}")
            return default
    
    return current


# =============================================================================
# 列过滤配置 - 定义需要保留的列和顺序
# =============================================================================

DESIRED_COLUMNS = [
    # 基础变异信息
    'Chr', 'Start', 'Ref', 'Alt', 
    # 基因功能注释
    'Func.refGene', 'Gene.refGene', 'HGNC_ID', 'HGVSc', 'Vcf_mut', 
    'GeneDetail.refGene', 'ExonicFunc.refGene', 'AAChange.refGene',
    # 数据库频率信息
    'avsnp151', 'EXON', 'INTRON', 'cytoBand', 'AF', 'AF_eas', 
    '1000g2015aug_all', '1000g2015aug_eas', 'esp6500siv2_all', 'ExAC_ALL',
    # 功能预测
    'SIFT_pred', 'LRT_pred', 'MetaSVM_pred', 'MetaLR_pred', 'MutationTaster_pred', 
    'Spliceai_masked',
    # 样本基因型信息
    'GT', 'AD', 'DP', 'GQ',
    # 临床注释
    'CLNSIG', 'CLNDN', 'CLNREVSTAT',
    # VEP注释结果
    'Consequence', 'Func.HGVS', 'MaxEntScan_Reduced_value'
]


# =============================================================================
# 配置初始化
# =============================================================================

# 默认配置（当没有提供配置文件时使用）
DEFAULT_CONFIG = {
    # ANNOVAR相关路径
    "annovar_script": "/haplox/users/chenjh/database/annovar/table_annovar.pl",
    "annovar_db": "/haplox/users/chenjh/database/annovar/humandb",
    
    # VEP相关路径
    "vep_env": "/haplox/users/chenjh/miniforge3/envs/vep_clone",
    "fasta_hs37d5": "/haplox/users/chenjh/database/genomes/hs37d5/hs37d5.chr.fa",
    "vep_cache_dir": "/haplox/users/chenjh/database/vep_cache",
    "vep_plugins_dir": "/haplox/users/chenjh/database/vep_plugins/",
    
    # 数据库文件路径
    "hgnc_transcript_db": "/haplox/users/chenjh/database/haima_file/hgnc_Transcript.txt",
    "exclude_chr_list": "/haplox/users/chenjh/database/haima_file/rm_chr.list",
    "spliceai_snv": "/haplox/users/chenjh/database/vep_plugins/spliceai/spliceai_scores.masked.snv.hg19.vcf.gz",
    "spliceai_indel": "/haplox/users/chenjh/database/vep_plugins/spliceai/spliceai_scores.masked.indel.hg19.vcf.gz",
    "maxentscan_dir": "/haplox/users/chenjh/database/vep_plugins/MaxEntScan",
    
    # 默认参数
    "ref": "hg19",
    "vep_cache_version": "115",  
    "vep_fork": "20"
}

# 全局配置变量，将在main函数中初始化
CONFIG = {}


# =============================================================================
# 工具函数模块
# =============================================================================

def get_Reduced_value_MaxEntScan(MaxEntScan_ref, MaxEntScan_alt):
    """
    计算MaxEntScan缩减值
    
    公式: (ref_score - alt_score) / ref_score
    
    参数:
        MaxEntScan_ref: 参考序列分数
        MaxEntScan_alt: 变异序列分数
        
    返回:
        float或str: 计算结果，无效输入返回'-'
    """
    if MaxEntScan_alt == '-' or MaxEntScan_ref == '-' or MaxEntScan_ref == 0:
        return '-'
    return (MaxEntScan_ref - MaxEntScan_alt) / MaxEntScan_ref


def str2float(val):
    """
    安全字符串转浮点数
    
    参数:
        val: 输入值
        
    返回:
        float或str: 转换后的浮点数，'-'保持不变
    """
    if val == '-':
        return '-'
    return float(val)


def convert_dataframe_to_json_serializable(df):
    """
    将DataFrame转换为JSON可序列化的格式
    
    处理步骤:
    1. 将NaN值转换为None
    2. 将numpy类型转换为Python原生类型
    3. 确保所有值都可以被JSON序列化
    
    参数:
        df: pandas DataFrame
        
    返回:
        list: 包含所有行的字典列表
    """
    # 复制DataFrame以避免修改原始数据
    df_copy = df.copy()
    
    # 将NaN值转换为None
    df_copy = df_copy.where(pd.notnull(df_copy), None)
    
    # 转换numpy类型为Python原生类型
    for col in df_copy.columns:
        if df_copy[col].dtype == np.int64:
            df_copy[col] = df_copy[col].astype(object).where(df_copy[col].notnull(), None)
        elif df_copy[col].dtype == np.float64:
            df_copy[col] = df_copy[col].astype(object).where(df_copy[col].notnull(), None)
    
    # 转换为字典列表
    records = df_copy.to_dict('records')
    
    return records


def filter_dataframe_columns(df, desired_columns):
    """
    过滤DataFrame，只保留指定的列并按照指定顺序排列
    
    参数:
        df: 原始DataFrame
        desired_columns: 需要保留的列名列表，按顺序排列
        
    返回:
        DataFrame: 过滤和重排序后的DataFrame
    """
    print("过滤和重排列数据...")
    
    # 检查哪些列在DataFrame中存在
    existing_columns = [col for col in desired_columns if col in df.columns]
    missing_columns = [col for col in desired_columns if col not in df.columns]
    
    if missing_columns:
        print(f"警告: 以下列在数据中不存在: {missing_columns}")
    
    print(f"保留 {len(existing_columns)} 列，按照指定顺序排列")
    
    # 只保留存在的列并按照指定顺序排列
    filtered_df = df[existing_columns].copy()
    
    return filtered_df


# =============================================================================
# VEP注释相关功能
# =============================================================================

def make_vep_input_from_annovar(input_path, vep_input_path):
    """
    将ANNOVAR多注释表转换为VEP输入格式
    
    处理步骤:
    1. 读取ANNOVAR输出文件
    2. 清理列名中的.1后缀
    3. 构建VEP所需的输入格式: Chr, start, End, ref/alt, ori, label
    4. 对于缺失(ref为'-')的情况调整起始位置
    
    参数:
        input_path: ANNOVAR输出文件路径
        vep_input_path: VEP输入文件输出路径
    """
    print("正在准备VEP输入文件...")
    
    # 分块读取ANNOVAR结果，避免内存问题
    df_annovar = pd.concat(
        pd.read_csv(input_path, header=0, sep="\t", chunksize=10000, low_memory=False),
        ignore_index=True
    )

    # 清理列名中的.1后缀
    df_annovar.columns = [c.replace(".1", "") for c in df_annovar.columns]

    # 提取VEP需要的列并格式化
    input_df_tovep = df_annovar[['Chr', 'Start', 'End', 'Ref', 'Alt']].copy()
    input_df_tovep['ori'] = '+'  # 默认正向链
    input_df_tovep['refalt'] = input_df_tovep.apply(lambda x: f"{x['Ref']}/{x['Alt']}", axis=1)
    
    # 对于缺失(ref为'-')的情况，起始位置+1
    input_df_tovep['start'] = input_df_tovep.apply(
        lambda x: int(x['Start']) + 1 if str(x['Ref']) == '-' else x['Start'], axis=1
    )
    
    input_df_tovep = input_df_tovep.drop(columns=['Ref', 'Alt'])
    input_df_tovep['row_number'] = input_df_tovep.index
    
    # 构建唯一标识符
    input_df_tovep['label'] = (
        input_df_tovep['Chr'].astype(str) + '_' +
        input_df_tovep['start'].astype(str) + '_' +
        input_df_tovep['End'].astype(str) + '_' +
        input_df_tovep['row_number'].astype(str)
    )
    
    # 选择VEP需要的列并输出
    input_df_tovep = input_df_tovep[['Chr', 'start', 'End', 'refalt', 'ori', 'label']]
    input_df_tovep.to_csv(vep_input_path, index=False, header=False, sep='\t')
    print(f"VEP输入文件已生成: {vep_input_path}")


def run_vep_annotation(input_file, vep_input, vep_output):
    """
    运行VEP注释
    
    处理步骤:
    1. 生成VEP输入文件
    2. 执行VEP注释命令
    
    参数:
        input_file: ANNOVAR输出文件路径
        vep_input: VEP输入文件路径
        vep_output: VEP输出文件路径
        
    返回:
        bool: 成功返回True，失败返回False
    """
    print("开始VEP注释...")

    # 生成VEP输入文件
    make_vep_input_from_annovar(input_file, vep_input)

    # 构建VEP命令
    cmd = [
        "vep",  # 直接使用系统环境中的vep命令
        "--fasta", CONFIG['fasta_hs37d5'],
        "--tab",
        "--fork", CONFIG['vep_fork'],
        "--force_overwrite",
        "--assembly", "GRCh37",
        "--merged",
        "--shift_hgvs=1",
        "--hgvs",
        "--use_given_ref",
        "--no_escape",
        "--offline",
        "--cache_version", CONFIG['vep_cache_version'],
        "--dir_cache", CONFIG['vep_cache_dir'],
        "--dir_plugins", CONFIG['vep_plugins_dir'],
        "--numbers",
        "--plugin", f"SpliceAI,snv={CONFIG['spliceai_snv']},indel={CONFIG['spliceai_indel']}",
        "--plugin", f"MaxEntScan,{CONFIG['maxentscan_dir']}",
        "-i", vep_input,
        "-o", vep_output,
        "--symbol",
        "--hgvs",
        "--canonical",
        "--flag_pick",
        "--pick_order", "refseq,canonical,appris,tsl,biotype,ccds,rank,length",
    ]

    # 执行VEP注释
    print(f"执行VEP命令: {' '.join(cmd)}")
    start = datetime.datetime.now()
    
    try:
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        stdout, stderr = proc.communicate()
        elapsed = (datetime.datetime.now() - start).total_seconds()

        if proc.returncode == 0:
            print(f"VEP执行成功！耗时 {elapsed:.1f}s")
            return True
        else:
            print("VEP执行失败：")
            if stdout:
                print(f"STDOUT: {stdout}")
            if stderr:
                print(f"STDERR: {stderr}")
            return False
    except Exception as e:
        print(f"VEP执行异常: {e}")
        return False


# =============================================================================
# ANNOVAR整合相关功能 - 拆分的子函数
# =============================================================================

def _load_transcript_database():
    """
    加载转录本数据库
    
    返回:
        DataFrame: 处理后的转录本数据库
    """
    print("加载转录本数据库...")
    hgncDB_df = pd.read_csv(CONFIG['hgnc_transcript_db'], sep='\t')
    
    # 重命名列以匹配脚本期望的列名
    hgncDB_df = hgncDB_df.rename(columns={
        'hgnc_id': 'HGNC_ID',
        'symbol': 'Approved_symbol', 
        'refseq_accession': 'MANE_Select_RefSeq'
    })
    
    # 处理HGNC_ID，去除"HGNC:"前缀
    hgncDB_df['HGNC_ID'] = hgncDB_df['HGNC_ID'].str.replace('HGNC:', '', regex=False)
    
    return hgncDB_df


def _process_vep_data(vep_df, hgncDB_df):
    """
    处理VEP注释数据
    
    参数:
        vep_df: VEP原始数据
        hgncDB_df: 转录本数据库
        
    返回:
        DataFrame: 处理后的VEP数据
    """
    print("处理VEP注释结果...")
    
    # 提取SpliceAI最大值（如果列存在）
    if 'SpliceAI_pred' in vep_df.columns:
        vep_df['Spliceai_masked'] = vep_df['SpliceAI_pred'].apply(
            lambda x: max([float(i) for i in x.split('|')[1:5]]) if x != '-' else '-'
        )
    else:
        vep_df['Spliceai_masked'] = '-'
        print("警告: SpliceAI_pred列不存在，Spliceai_masked设为默认值'-'")
    
    # 提取ANNOVAR顺序
    vep_df['annovar_order'] = vep_df['#Uploaded_variation'].apply(lambda x: int(x.split('_')[3]))
    
    # 去除转录本版本号
    if 'Feature' in vep_df.columns:
        vep_df['NM'] = vep_df['Feature'].apply(lambda x: x.split('.')[0] if x != '-' else '-')
    else:
        vep_df['NM'] = '-'
        print("警告: Feature列不存在，NM设为默认值'-'")
    
    # 构建基因ID和名称映射
    valid_id_hgncDB_df = hgncDB_df[hgncDB_df['HGNC_ID'] != '-']
    hgncDB_ts2id = dict(zip(valid_id_hgncDB_df['MANE_Select_RefSeq'], valid_id_hgncDB_df['HGNC_ID']))
    hgncDB_id2gene = dict(zip(valid_id_hgncDB_df['HGNC_ID'], valid_id_hgncDB_df['Approved_symbol']))
    
    # 修正VEP中的基因ID和名称
    print("修正基因注释信息...")
    if 'HGNC_ID' in vep_df.columns:
        vep_df.loc[vep_df['HGNC_ID'] == '-', 'HGNC_ID'] = vep_df[vep_df['HGNC_ID'] == '-']['NM'].apply(
            lambda x: hgncDB_ts2id.get(x, '-')
        )
    else:
        vep_df['HGNC_ID'] = vep_df['NM'].apply(lambda x: hgncDB_ts2id.get(x, '-'))
    
    if 'SYMBOL' in vep_df.columns:
        vep_df.loc[vep_df['HGNC_ID'] != '-', 'SYMBOL'] = vep_df[vep_df['HGNC_ID'] != '-']['HGNC_ID'].apply(
            lambda x: hgncDB_id2gene.get(str(x), '-')
        )
    else:
        vep_df['SYMBOL'] = vep_df['HGNC_ID'].apply(lambda x: hgncDB_id2gene.get(str(x), '-'))
    
    return vep_df, hgncDB_ts2id, hgncDB_id2gene


def _select_best_transcript(vep_df, hgncDB_df):
    """
    选择最优转录本
    
    优先级: MANE > PICK > 第一个
    
    参数:
        vep_df: 处理后的VEP数据
        hgncDB_df: 转录本数据库
        
    返回:
        DataFrame: 选择最优转录本后的数据
    """
    print("选择最优转录本...")
    
    # 标记MANE转录本
    if 'PICK' in vep_df.columns:
        vep_df.loc[vep_df['NM'].isin(hgncDB_df['MANE_Select_RefSeq']), 'PICK'] = '2'
    
    # 按优先级选择转录本：MANE > PICK > 第一个
    vep_df_1_mane = vep_df[vep_df['PICK'] == '2'] if 'PICK' in vep_df.columns else pd.DataFrame()
    mane_id = vep_df_1_mane['#Uploaded_variation'] if not vep_df_1_mane.empty else []
    filtered_df = vep_df[~vep_df['#Uploaded_variation'].isin(mane_id)]
    
    vep_df_2_pick = filtered_df[filtered_df['PICK'] == '1'] if 'PICK' in vep_df.columns else pd.DataFrame()
    pick_id = vep_df_2_pick['#Uploaded_variation'] if not vep_df_2_pick.empty else []
    filtered_df = filtered_df[~filtered_df['#Uploaded_variation'].isin(pick_id)]
    
    unique_df = filtered_df.drop_duplicates(subset='#Uploaded_variation')  # 第一个转录本
    
    # 合并所有选择的转录本
    combined_df = pd.concat([vep_df_1_mane, vep_df_2_pick, unique_df], ignore_index=True)
    sorted_df = combined_df.sort_values(by='annovar_order')
    vep_df = sorted_df.reset_index(drop=True)
    
    return vep_df


def _process_annovar_data(annovar_df):
    """
    处理ANNOVAR注释数据
    
    参数:
        annovar_df: ANNOVAR原始数据
        
    返回:
        DataFrame: 处理后的ANNOVAR数据
    """
    print("处理ANNOVAR注释结果...")
    
    # 构建VCF突变格式
    annovar_df['Vcf_mut'] = (annovar_df['Otherinfo4'].astype(str) + ':' + 
                           annovar_df['Otherinfo5'].astype(str) + ':' + 
                           annovar_df['Otherinfo7'].astype(str) + '/' + 
                           annovar_df['Otherinfo8'].astype(str))
    
    # 提取基因型信息
    if 'Otherinfo12' in annovar_df.columns and len(annovar_df) > 0:
        if annovar_df['Otherinfo12'].iloc[0].split(':')[0:4] == ['GT', 'AD', 'DP', 'GQ']:
            genotype_cols = annovar_df['Otherinfo12'].iloc[0].split(':')[0:4]
            annovar_df[genotype_cols] = annovar_df['Otherinfo13'].str.split(':', expand=True).iloc[:, 0:4]
            drop_col = ["Otherinfo1","Otherinfo2","Otherinfo3","Otherinfo4","Otherinfo5","Otherinfo6",
                      "Otherinfo7","Otherinfo8","Otherinfo9","Otherinfo10","Otherinfo11","Otherinfo12","Otherinfo13"]
        elif annovar_df['Otherinfo12'].iloc[0].split(':')[0:4] == ['GT', 'AD', 'AF', 'DP']:
            genotype_cols = annovar_df['Otherinfo12'].iloc[0].split(':')[0:4]
            annovar_df[genotype_cols] = annovar_df['Otherinfo14'].str.split(':', expand=True).iloc[:, 0:4]
            annovar_df['GQ'] = ""
            drop_col = ["Otherinfo1","Otherinfo2","Otherinfo3","Otherinfo4","Otherinfo5","Otherinfo6",
                      "Otherinfo7","Otherinfo8","Otherinfo9","Otherinfo10","Otherinfo11","Otherinfo12",
                      "Otherinfo13","Otherinfo14"]
        else:
            drop_col = [col for col in annovar_df.columns if col.startswith('Otherinfo')]
    else:
        drop_col = [col for col in annovar_df.columns if col.startswith('Otherinfo')]
    
    annovar_df = annovar_df.drop(drop_col, axis=1)
    
    # 标记InDel类型
    annovar_df.loc[(annovar_df['Ref'] == '-') | (annovar_df['Alt'] == '-'), 'Mut_type'] = 'InDel'
    
    # 格式化基因型信息
    if 'GT' in annovar_df.columns:
        annovar_df['GT'] = ' ' + annovar_df['GT'].astype(str)
    if 'AD' in annovar_df.columns:
        annovar_df['AD'] = ' ' + annovar_df['AD'].str.replace(",", "/")
    
    return annovar_df


def _merge_annotations(annovar_df, vep_df, hgncDB_id2gene, hgncDB_df_copy):
    """
    合并ANNOVAR和VEP注释结果
    
    参数:
        annovar_df: 处理后的ANNOVAR数据
        vep_df: 处理后的VEP数据
        hgncDB_id2gene: 基因ID到名称的映射
        hgncDB_df_copy: 转录本数据库副本
        
    返回:
        DataFrame: 合并后的注释结果
    """
    print("合并VEP和ANNOVAR注释...")
    
    # 定义要添加的列
    add_dict = ['Spliceai_masked']
    
    # 添加可能存在的VEP列
    optional_columns = ["MaxEntScan_alt", "MaxEntScan_diff", "MaxEntScan_ref",
                       "Consequence", "SYMBOL", "EXON", "INTRON", "HGVSc", "HGVSp", "HGNC_ID"]
    
    for col in optional_columns:
        if col in vep_df.columns:
            add_dict.append(col)
    
    # 添加空列
    for key in add_dict:
        annovar_df[key] = ""
    
    # 按顺序合并数据
    annovar_df_new = annovar_df.iloc[list(vep_df['annovar_order'].values)]
    annovar_df_new = annovar_df_new.reset_index(drop=True)
    annovar_df_new.loc[:, add_dict] = vep_df.loc[:, add_dict]
    
    # ==================== 额外标注和格式化 ====================
    print("添加额外注释信息...")
    
    # 添加HGNC基因名称
    annovar_df_new['HGNC'] = annovar_df_new['HGNC_ID'].apply(lambda x: hgncDB_id2gene.get(x, '-'))
    annovar_df_new['Func.HGVS'] = annovar_df_new['Consequence'] if 'Consequence' in annovar_df_new.columns else ''
    
    # 处理缺失值
    if 'HGVSc' in annovar_df_new.columns:
        annovar_df_new.loc[annovar_df_new['HGVSc'] == "-", 'HGVSc'] = 'NA'
    if 'HGVSp' in annovar_df_new.columns:
        annovar_df_new.loc[annovar_df_new['HGVSp'] == "-", 'HGVSp'] = 'NA'
    
    # 处理MaxEntScan数值（如果列存在）
    for _col in ['MaxEntScan_ref', 'MaxEntScan_diff', 'MaxEntScan_alt']:
        if _col in annovar_df_new.columns:
            annovar_df_new[_col] = annovar_df_new[_col].map(str2float)
    
    # 计算MaxEntScan缩减值（如果相关列存在）
    if all(col in annovar_df_new.columns for col in ['MaxEntScan_ref', 'MaxEntScan_alt']):
        annovar_df_new['MaxEntScan_Reduced_value'] = annovar_df_new.apply(
            lambda x: get_Reduced_value_MaxEntScan(x['MaxEntScan_ref'], x['MaxEntScan_alt']), axis=1
        )
    else:
        annovar_df_new['MaxEntScan_Reduced_value'] = '-'
    
    # ==================== 过滤和最终处理 ====================
    print("过滤和转换最终结果...")
    
    # 只保留SYMBOL与Gene.refGene匹配的行
    if 'SYMBOL' in annovar_df_new.columns and 'Gene.refGene' in annovar_df_new.columns:
        annovar_df_new = annovar_df_new[annovar_df_new['SYMBOL'] == annovar_df_new['Gene.refGene']]
    
    # 转换HGNC_ID为NM编号 - 使用完整的MANE_Select_RefSeq
    hgnc_dict_NM = dict(zip(hgncDB_df_copy['HGNC_ID'].astype(str), hgncDB_df_copy['MANE_Select_RefSeq']))
    annovar_df_new['HGNC_ID'] = annovar_df_new['HGNC_ID'].apply(lambda x: hgnc_dict_NM.get(x, '-'))
    
    return annovar_df_new


# =============================================================================
# 主整合函数
# =============================================================================

def run_annovar_integration(vep_file, annovar_file, sample, outdir, filter_columns=True):
    """
    整合ANNOVAR和VEP注释结果
    
    参数:
        vep_file: VEP输出文件路径
        annovar_file: ANNOVAR输出文件路径
        sample: 样本名称
        outdir: 输出目录
        filter_columns: 是否进行列筛选，默认为True
        
    返回:
        bool: 成功返回True，失败返回False
    """
    try:
        # 读取输入文件
        print("读取VEP和ANNOVAR结果文件...")
        
        # 首先检查VEP文件的实际列名
        print("检查VEP文件列名...")
        with open(vep_file, 'r') as f:
            for i, line in enumerate(f):
                if i == 59:  # 找到列名行
                    actual_columns = line.strip().split('\t')
                    print(f"VEP文件实际列名: {actual_columns}")
                    break
        
        # 动态构建需要的列名列表
        usecols = []
        expected_columns = [
            "#Uploaded_variation", "Consequence", "SYMBOL", "Feature", "HGNC_ID", "PICK",
            "EXON", "INTRON", "HGVSc", "HGVSp", "HGVS_OFFSET", "SpliceAI_pred",
            "MaxEntScan_alt", "MaxEntScan_diff", "MaxEntScan_ref"
        ]
        
        # 检查哪些列实际存在
        available_columns = []
        for col in expected_columns:
            if col in actual_columns:
                available_columns.append(col)
            else:
                print(f"警告: 列 '{col}' 在VEP输出中不存在")
        
        print(f"可用的列: {available_columns}")
        
        # 读取VEP文件，只使用存在的列
        vep_df = pd.read_csv(vep_file, sep='\t', low_memory=False, header=59, usecols=available_columns)
        
        # 读取ANNOVAR文件
        annovar_df = pd.read_csv(annovar_file, sep='\t', low_memory=False)
        
        # 加载转录本数据库
        hgncDB_df = _load_transcript_database()
        
        # 创建输出目录
        os.makedirs(outdir, exist_ok=True)
        final_output = f"{outdir}/{sample}_anno.txt"
        json_output = f"{outdir}/{sample}_anno.json"
        
        # 处理转录本库
        hgncDB_df['HGNC_ID'] = hgncDB_df['HGNC_ID'].astype(str)
        hgncDB_df_copy = hgncDB_df.copy()
        
        # 处理VEP数据
        vep_df, hgncDB_ts2id, hgncDB_id2gene = _process_vep_data(vep_df, hgncDB_df)
        
        # 选择最优转录本
        vep_df = _select_best_transcript(vep_df, hgncDB_df)
        
        # 处理ANNOVAR数据
        annovar_df = _process_annovar_data(annovar_df)
        
        # 合并注释结果
        annovar_df_new = _merge_annotations(annovar_df, vep_df, hgncDB_id2gene, hgncDB_df_copy)
        
        # ==================== 输出结果 ====================
        print("输出最终结果...")
        
        # 输出完整的TXT格式
        annovar_df_new.to_csv(final_output, index=None, sep='\t')
        print(f"完整TXT格式结果: {final_output}")
        
        # 输出JSON格式
        print("生成JSON格式结果...")
        
        # 根据filter_columns参数决定是否进行列筛选
        if filter_columns:
            print("进行列筛选...")
            # 过滤列并重新排序
            output_df = filter_dataframe_columns(annovar_df_new, DESIRED_COLUMNS)
        else:
            print("跳过列筛选，保留所有列...")
            output_df = annovar_df_new.copy()
        
        # 准备JSON数据结构
        json_data = {
            "sample_id": sample,
            "pipeline_version": "1.0",
            "timestamp": datetime.datetime.now().isoformat(),
            "total_variants": len(output_df),
            "filter_columns": filter_columns,
            "columns": list(output_df.columns),
            "data": convert_dataframe_to_json_serializable(output_df)
        }
        
        # 写入JSON文件
        with open(json_output, 'w', encoding='utf-8') as f:
            json.dump(json_data, f, indent=2, ensure_ascii=False)
        
        print(f"JSON格式结果: {json_output}")
        print(f"ANNOVAR和VEP结果整合完成")
        return True
        
    except Exception as e:
        print(f"整合过程中出现错误: {e}")
        import traceback
        traceback.print_exc()
        return False


# =============================================================================
# 主流程功能
# =============================================================================

def run_annovar_annotation(input_vcf, output_prefix, ref="hg19"):
    """
    运行ANNOVAR注释
    
    参数:
        input_vcf: 输入VCF文件路径
        output_prefix: 输出文件前缀
        ref: 基因组版本
        
    返回:
        bool: 成功返回True，失败返回False
    """
    print("开始ANNOVAR注释...")
    
    cmd = [
        "perl", CONFIG['annovar_script'],
        input_vcf,
        CONFIG['annovar_db'],
        "-buildver", ref,
        "-out", output_prefix,
        "-remove",
        "-polish",
        "-protocol", "refGene,cytoBand,esp6500siv2_all,exac03,gnomad211_genome,avsnp151,dbnsfp42c,clinvar_20250721,1000g2015aug_all,1000g2015aug_eas",
        "-operation", "g,r,f,f,f,f,f,f,f,f",
        "-arg", "-splicing_threshold 10 -hgvs,,,,,,,,,",
        "-nastring", ".",
        "-vcfinput"
    ]    
    print(f"执行ANNOVAR命令: {' '.join(cmd)}")
    
    try:
        result = subprocess.run(cmd, check=True, capture_output=True, text=True)
        print("ANNOVAR注释完成")
        return True
    except subprocess.CalledProcessError as e:
        print(f"ANNOVAR注释失败: {e}")
        if e.stderr:
            print(f"STDERR: {e.stderr}")
        return False
    except Exception as e:
        print(f"ANNOVAR执行异常: {e}")
        return False


def filter_chromosomes(input_file, output_file):
    """
    过滤掉不需要的染色体和MT变异位点
    
    参数:
        input_file: 输入文件路径
        output_file: 输出文件路径
        
    返回:
        bool: 成功返回True，失败返回False
    """
    print("过滤染色体...")
    
    try:
        # 读取要排除的染色体列表
        with open(CONFIG['exclude_chr_list'], 'r') as f:
            exclude_chrs = set(line.strip() for line in f)
        
        # 读取输入文件
        df = pd.read_csv(input_file, sep='\t', low_memory=False)
        
        # 过滤条件：不在排除列表中且不是chrMT
        mask = (~df['Chr'].astype(str).isin(exclude_chrs)) & (~df['Chr'].astype(str).str.upper().str.startswith('CHRMT'))
        
        # 应用过滤
        filtered_df = df[mask]
        
        # 保存结果
        filtered_df.to_csv(output_file, sep='\t', index=False)
        print(f"染色体过滤完成: {len(filtered_df)} 行保留 (原始: {len(df)} 行)")
        return True
        
    except Exception as e:
        print(f"染色体过滤失败: {e}")
        return False


# =============================================================================
# 主函数
# =============================================================================

def main():
    """
    主函数 - 协调整个注释流程
    
    流程概述:
        1. ANNOVAR注释
        2. 染色体过滤  
        3. VEP注释
        4. 结果整合和输出
    """
    parser = argparse.ArgumentParser(
        description='整合的ANNOVAR和VEP注释流程',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
使用示例:
  python unified_annotation_pipeline.py -vcf input.vcf -s sample001 -o ./results
  python unified_annotation_pipeline.py --vcf input.vcf --sample sample001 --output ./results --ref hg38
        """
    )
    
    parser.add_argument('-vcf', '--vcf', required=True, help='输入VCF文件路径')
    parser.add_argument('-s', '--sample', required=True, help='样本名称')
    parser.add_argument('-o', '--output', required=True, help='输出目录路径')
    parser.add_argument('--ref', default=None, 
                       help='基因组版本 (默认: 从配置文件读取)')
    parser.add_argument('--annovar_db', default=None,
                       help='ANNOVAR数据库路径 (默认: 从配置文件读取)')
    # 添加控制列筛选的参数
    parser.add_argument('--no_filter_columns', action='store_true', 
                       help='不进行列筛选，保留所有列 (默认进行列筛选)')
    # 添加配置文件参数
    parser.add_argument('--config', required=True, help='YAML配置文件路径')
    
    args = parser.parse_args()
    
    # 验证输入文件
    if not os.path.exists(args.vcf):
        print(f"错误: VCF文件不存在: {args.vcf}")
        sys.exit(1)
    
    # 加载配置文件
    yaml_config = load_config(args.config)
    
    # 初始化全局配置
    global CONFIG
    CONFIG = {}
    
    # 从YAML配置中获取annotation部分的配置
    annotation_config = get_config_value(yaml_config, 'annotation', {})
    
    # 更新CONFIG：先使用默认配置，然后用YAML配置覆盖
    CONFIG.update(DEFAULT_CONFIG)
    CONFIG.update(annotation_config)
    
    # 处理命令行参数覆盖配置的情况
    if args.ref is not None:
        CONFIG['ref'] = args.ref
    if args.annovar_db is not None:
        CONFIG['annovar_db'] = args.annovar_db
    
    # 创建输出目录
    os.makedirs(args.output, exist_ok=True)
    
    print("=" * 60)
    print("ANNOVAR和VEP整合注释流程")
    print("=" * 60)
    print(f"样本: {args.sample}")
    print(f"输入VCF: {args.vcf}")
    print(f"输出目录: {args.output}")
    print(f"基因组版本: {CONFIG['ref']}")
    print(f"配置文件: {args.config}")
    print(f"列筛选: {'否' if args.no_filter_columns else '是'}")
    print(f"开始时间: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("-" * 60)
    
    # 步骤1: 运行ANNOVAR注释
    annovar_prefix = os.path.join(args.output, f"{args.sample}.annovar")
    if not run_annovar_annotation(args.vcf, annovar_prefix, CONFIG['ref']):
        print("ANNOVAR注释失败，流程终止")
        sys.exit(1)
    
    # ANNOVAR输出文件
    annovar_output = f"{annovar_prefix}.{CONFIG['ref']}_multianno.txt"
    final_annovar_output = os.path.join(args.output, f"{args.sample}.filter.{CONFIG['ref']}_multianno.txt")
    
    # 步骤2: 过滤染色体
    if not filter_chromosomes(annovar_output, final_annovar_output):
        print("染色体过滤失败，流程终止")
        sys.exit(1)
    
    # 步骤3: 运行VEP注释
    vep_input = os.path.join(args.output, f"{args.sample}_vep_input")
    vep_output = os.path.join(args.output, f"{args.sample}_vep_output")
    
    if not run_vep_annotation(final_annovar_output, vep_input, vep_output):
        print("VEP注释失败，流程终止")
        sys.exit(1)
    
    # 清理VEP临时文件
    vep_summary = f"{vep_output}_summary.html"
    for temp_file in [vep_input, vep_summary]:
        if os.path.exists(temp_file):
            os.remove(temp_file)
            print(f"已删除临时文件: {temp_file}")
    
    # 步骤4: 整合ANNOVAR和VEP结果
    # 根据参数决定是否进行列筛选
    filter_columns = not args.no_filter_columns
    
    if not run_annovar_integration(vep_output, final_annovar_output, args.sample, args.output, filter_columns):
        print("结果整合失败，流程终止")
        sys.exit(1)
    
    # 最终输出文件
    final_output = os.path.join(args.output, f"{args.sample}_anno.txt")
    json_output = os.path.join(args.output, f"{args.sample}_anno.json")
    
    print("=" * 60)
    print("流程完成!")
    print(f"完整注释结果 (TXT): {final_output}")
    print(f"注释结果 (JSON): {json_output}")
    print(f"列筛选: {'是' if filter_columns else '否'}")
    print(f"完成时间: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("=" * 60)


if __name__ == "__main__":
    main()