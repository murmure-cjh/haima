#!/usr/bin/env python3
# encoding: utf-8
"""
外显子测序数据质量控制分析脚本

创建时间: 2025.11.25
最后修改时间: 2025.11.25
版本: v1.0.0
作者: chenjh

描述：
    外显子测序数据质量控制分析脚本，用于统计处理软件fastp和bamdst软件的结果，
    生成JSON和TXT格式的质控报告，以及Shell脚本要求的CSV报告。

主要功能：
    1. 解析dedup_metrics文件获取基本测序指标
    2. 解析panel区域的coverage报告
    3. 解析panel区域的深度分布并计算FLOAD指标
    4. 解析fastp JSON文件获取质量指标
    5. 生成JSON和TXT格式的质控报告
    6. 生成Shell脚本要求的CSV报告

输入参数：
    sample: 样本名称
    bamfile: BAM文件路径
    dedup_metrics: 重复标记质控结果文件路径
    panelfile: panel区域文件路径
    outdir: 输出目录
    sampleid: 样本ID

输出文件：
    {sample}.sample_information.json: JSON格式质控报告
    {sample}.sample_information.txt: TXT格式质控报告
    {sample}_QC_result.csv: Shell脚本要求的CSV报告

基础来源：
    /x05_haplox/users/wangx/xiaohaima/necessary_files/database/QC_fload_V5.pl
    /x05_haplox/users/wangx/hereditary/miniwes_pipeline/QC_result/parse_QC_result.sh
"""

import os
import json
import sys
import csv
from pathlib import Path
from typing import Dict, Any, List



def qc_fload_analysis(sample, bamfile, dedup_metrics, panelfile, outdir, sampleid):
    """
    外显子测序数据质量控制分析主函数
    
    参数:
        sample: 样本名称
        bamfile: BAM文件路径
        dedup_metrics: 重复标记质控结果文件路径
        panelfile: panel区域文件路径
        outdir: 输出目录路径
        sampleid: 样本ID
    
    返回:
        None
    """
    
    # 初始化所有质控变量
    qc_variables = initialize_qc_variables()
    
    # 创建临时目录用于存储中间文件
    temp_directory = Path(outdir) / "tmp"
    temp_directory.mkdir(exist_ok=True)
    
    # 创建result目录用于存储Shell脚本要求的报告
    result_directory = Path(outdir) 
    result_directory.mkdir(exist_ok=True)
    
    try:
        # ==================== 数据解析阶段 ====================
        
        # 1. 解析dedup_metrics文件获取基本测序指标
        parse_dedup_metrics_file(dedup_metrics, sample, qc_variables)
        print(Path(outdir))
        # 2. 解析panel区域的coverage报告
        parse_coverage_report_file(Path(outdir) / "coverage.report", "panel", qc_variables)
        
        # 3. 解析panel区域的深度分布并计算FLOAD指标
        parse_depth_distribution_file(Path(outdir) / "depth_distribution.plot", "panel", qc_variables)
        
        # 4. 解析fastp JSON文件获取测序质量指标
        fastp_quality_data = parse_fastp_json_file(Path(outdir) / f"{sample}_cfDNA_fastp.json")
        
        # ==================== 报告生成阶段 ====================
        
        # 5. 生成JSON格式的质控报告
        generate_json_qc_report(outdir, sample, qc_variables, fastp_quality_data, sampleid)
        
        # 6. 生成TXT格式的质控报告
        generate_txt_qc_report(outdir, sample, qc_variables, fastp_quality_data)
        
        # 7. 生成Shell脚本要求的CSV报告
        generate_shell_script_csv_report(Path(outdir), sample, qc_variables, fastp_quality_data)
                
        print(f"质控分析完成: {sample}")
        
    except Exception as error:
        print(f"质控分析失败: {error}")
        sys.exit(1)


def generate_shell_script_csv_report(result_dir: Path, sample: str, 
                                    qc_variables: Dict[str, Any], 
                                    fastp_quality_data: Dict[str, float]) -> None:
    """
    生成Shell脚本要求的CSV报告
    
    参数:
        result_dir: 结果目录路径
        sample: 样本名称
        qc_variables: 质控变量字典
        fastp_quality_data: fastp质量数据
    
    返回:
        None
    """
    # 生成CSV报告文件
    csv_file_path = result_dir / f"{sample}_QC_result.csv"
    
    # Shell脚本要求的表头
    header = "样本编号,目标区域覆盖度(%),平均测序深度(X),20X 覆盖度(%),30X 覆盖度(%),Q20(质量值≥20的碱基所占百分比)(%),Q30(质量值≥30的碱基所占百分比)(%)"
    
    # 提取Shell脚本需要的指标
    # 注意：Shell脚本中提取的是panel区域的数据
    mapping_rate = qc_variables['Mapping_rate']
    average_depth = qc_variables['panelavedepth']
    coverage_20x = qc_variables['cov20_rate']
    coverage_30x = qc_variables['cov30_rate']
    q20_rate = fastp_quality_data['q20_rate']
    q30_rate = fastp_quality_data['q30_rate']
    
    # 写入CSV文件
    with open(csv_file_path, 'w', encoding='utf-8-sig') as csv_file:  # utf-8-sig for Excel compatibility
        csv_file.write(header + '\n')
        csv_file.write(f"{sample},{mapping_rate},{average_depth},{coverage_20x},{coverage_30x},{q20_rate},{q30_rate}\n")
    
    print(f"Shell脚本CSV报告已生成: {csv_file_path}")


def initialize_qc_variables() -> Dict[str, Any]:
    """
    初始化所有质控变量
    
    返回:
        Dict[str, Any]: 包含所有质控变量的字典，初始值为0
    """
    return {
        # 基本测序指标
        'Toatl_clean_read': 0,      # 总clean reads数
        'Mapping_rate': 0,          # 比对率(%)
        'Dup_rate': 0,              # 重复率(%)
        
        # Panel区域指标
        'panelCapture_rate': 0,     # panel区域捕获率
        'panelavedepth': 0,         # panel区域平均深度
        'panelCoverage': 0,         # panel区域1X覆盖率
        'panelmap_read': 0,         # panel区域比对reads数(百万)
        'panelmap_base': 0,         # panel区域比对bases数(GB)
        'paneltarget_base': 0,      # panel区域目标bases数(GB)
        'panelcov10_rate': 0,       # panel区域10X覆盖率
        'panelcov100_rate': 0,      # panel区域100X覆盖率
        'fload80': 0,               # panel区域FLOAD80指标
        
        # 深度覆盖率指标
        'cov0_rate': 0,             # 0X覆盖率
        'cov5_rate': 0,             # 5X覆盖率
        'cov10_rate': 0,            # 10X覆盖率
        'cov20_rate': 0,            # 20X覆盖率
        'cov30_rate': 0,            # 30X覆盖率
        'cov50_rate': 0,            # 50X覆盖率
        'cov100_rate': 0,           # 100X覆盖率
    }


def parse_dedup_metrics_file(dedup_metrics_path: str, sample_name: str, 
                          qc_variables: Dict[str, Any]) -> None:
    """
    解析dedup_metrics文件获取基本测序指标
    
    参数:
        dedup_metrics_path: 重复标记质控结果文件路径
        sample_name: 样本名称
        qc_variables: 质控变量字典
    
    返回:
        None
    """
    try:
        with open(dedup_metrics_path, 'r', encoding='utf-8') as file:
            for line in file:
                line = line.strip()
                if line.startswith(sample_name):
                    fields = line.split()
                    
                    # 计算重复率
                    duplication_rate = float(fields[8]) * 100
                    qc_variables['Dup_rate'] = round(duplication_rate, 2)
                    
                    # 计算总clean reads数
                    total_clean_reads = (
                        int(fields[1]) + int(fields[3]) + 
                        (int(fields[2]) * 2) + int(fields[4])
                    )
                    qc_variables['Toatl_clean_read'] = total_clean_reads
                    
                    # 计算比对率
                    mapped_reads = total_clean_reads - int(fields[4])
                    mapping_rate = mapped_reads * 100 / total_clean_reads
                    qc_variables['Mapping_rate'] = round(mapping_rate, 2)
                    break
                    
    except Exception as error:
        print(f"解析dedup_metrics文件错误 {dedup_metrics_path}: {error}")


def parse_coverage_report_file(report_file_path, region_type, qc_variables):
    """
    解析coverage.report文件，提取Panel区域特定深度的覆盖度信息。
    Args:
        report_file_path: 报告文件路径
        region_type: 区域类型 (panel)
        qc_variables: 存储质控数据的字典
    """
    # 检查文件是否存在
    if not os.path.exists(str(report_file_path)):
        print(f"Warning: 文件不存在 {report_file_path}")
        return

    try:
        with open(report_file_path, 'r') as f:
            for line in f:
                line = line.strip()
                
                # 优化：只处理以 [Target] 开头的行
                if not line.startswith('[Target]'):
                    continue
                
                # 解析Average depth字段
                if "[Target] Average depth" in line:
                    try:
                        parts = line.split('\t')
                        if len(parts) >= 2:
                            value_str = parts[-1].strip()
                            qc_variables['panelavedepth'] = float(value_str)
                    except ValueError:
                        print(f"无法解析Average depth值: {line}")
                        continue
                
                # 解析Fraction of Target Data in mapped data字段
                elif "[Target] Fraction of Target Data in mapped data" in line:
                    try:
                        parts = line.split('\t')
                        if len(parts) >= 2:
                            value_str = parts[-1].strip()
                            # 处理可能存在的百分号
                            if value_str.endswith('%'):
                                value_str = value_str.rstrip('%')
                            # 直接使用原始值，不进行百分比转换
                            qc_variables['panelCapture_rate'] = float(value_str)
                    except ValueError:
                        print(f"无法解析Fraction of Target Data值: {line}")
                        continue
                
                # 解析覆盖度字段
                elif "Fraction Region covered" in line:
                    try:
                        parts = line.split('\t')
                        if len(parts) < 2: 
                            continue
                        
                        value_str = parts[-1].strip()
                        if value_str.endswith('%'):
                            value = float(value_str.strip('%')) / 100.0
                        else:
                            value = float(value_str)
                    except ValueError:
                        continue

                    # 赋值给panel区域变量
                    if "Fraction Region covered > 0x" in line or "Fraction Region covered (>=1x)" in line:
                        qc_variables['panelCoverage'] = value * 100  # 转换为百分比
                    elif "Fraction Region covered (>=5x)" in line:
                        qc_variables['cov5_rate'] = value * 100  # 转换为百分比
                    elif "Fraction Region covered (>=10x)" in line:
                        qc_variables['panelcov10_rate'] = value * 100  # 转换为百分比
                    elif "Fraction Region covered (>=20x)" in line:
                        qc_variables['cov20_rate'] = value * 100  # 转换为百分比
                    elif "Fraction Region covered (>=30x)" in line:
                        qc_variables['cov30_rate'] = value * 100  # 转换为百分比
                    elif "Fraction Region covered (>=50x)" in line:
                        qc_variables['cov50_rate'] = value * 100  # 转换为百分比
                    elif "Fraction Region covered (>=100x)" in line:
                        qc_variables['panelcov100_rate'] = value * 100  # 转换为百分比

    except Exception as error:
        print(f"解析coverage报告文件错误 {report_file_path}: {error}")


def parse_depth_distribution_file(plot_file_path: Path, region_type: str, 
                                 qc_variables: Dict[str, Any]) -> None:
    """
    解析深度分布文件并计算FLOAD指标
    
    参数:
        plot_file_path: 深度分布文件路径
        region_type: 区域类型 ('panel')
        qc_variables: 质控变量字典
    
    返回:
        None
    """
    if not plot_file_path.exists():
        print(f"警告: 文件不存在 {plot_file_path}")
        return
        
    try:
        total_depth_sum = 0
        depth_distribution_data = []
        
        # 第一遍读取：计算总深度和不同深度的覆盖率
        with open(plot_file_path, 'r', encoding='utf-8') as file:
            for line in file:
                line = line.strip()
                if not line:
                    continue
                    
                parts = line.split('\t')
                if len(parts) >= 5:
                    depth = int(parts[0])
                    depth_count = int(parts[1])
                    coverage_rate = float(parts[4])
                    
                    # 累加总深度
                    total_depth_sum += depth * depth_count
                    depth_distribution_data.append((depth, depth_count, coverage_rate))
                    
                    # 记录特定深度的覆盖率
                    update_coverage_rates(region_type, depth, coverage_rate, qc_variables)
        
        # 计算FLOAD 80指标
        calculate_fload_80(region_type, total_depth_sum, depth_distribution_data, qc_variables)
                
    except Exception as error:
        print(f"解析深度分布文件错误 {plot_file_path}: {error}")


def update_coverage_rates(region_type: str, depth: int, coverage_rate: float, 
                         qc_variables: Dict[str, Any]) -> None:
    """
    更新不同深度的覆盖率指标
    
    参数:
        region_type: 区域类型 ('panel')
        depth: 深度值
        coverage_rate: 覆盖率
        qc_variables: 质控变量字典
    
    返回:
        None
    """
    coverage_percentage = coverage_rate * 100
    
    if region_type == "panel":
        if depth == 0:
            qc_variables['cov0_rate'] = coverage_percentage
        elif depth == 5:
            qc_variables['cov5_rate'] = coverage_percentage
        elif depth == 10:
            qc_variables['panelcov10_rate'] = coverage_percentage
        elif depth == 20:
            qc_variables['cov20_rate'] = coverage_percentage
        elif depth == 30:
            qc_variables['cov30_rate'] = coverage_percentage
        elif depth == 50:
            qc_variables['cov50_rate'] = coverage_percentage
        elif depth == 100:
            qc_variables['panelcov100_rate'] = coverage_percentage


def calculate_fload_80(region_type: str, total_depth_sum: int, 
                      depth_distribution_data: List[tuple], 
                      qc_variables: Dict[str, Any]) -> None:
    """
    计算FLOAD 80指标
    
    参数:
        region_type: 区域类型 ('panel')
        total_depth_sum: 总深度和
        depth_distribution_data: 深度分布数据
        qc_variables: 质控变量字典
    
    返回:
        None
    """
    cumulative_depth_sum = 0
    fload80_depth = 0
    target_sum = total_depth_sum * 0.2  # FLOAD80对应的深度阈值
    
    # 计算累积深度，找到FLOAD80对应的深度
    for depth, depth_count, _ in depth_distribution_data:
        cumulative_depth_sum += depth * depth_count
        if cumulative_depth_sum > target_sum:
            fload80_depth = depth
            break
    
    # 计算FLOAD80值
    if fload80_depth > 0:
        if region_type == "panel":
            average_depth = qc_variables['panelavedepth']
            qc_variables['fload80'] = average_depth / fload80_depth


def parse_fastp_json_file(json_file_path: Path) -> Dict[str, float]:
    """
    解析fastp JSON文件获取质量指标
    
    参数:
        json_file_path: fastp JSON文件路径
    
    返回:
        Dict[str, float]: 包含fastp质量指标的字典
    """
    fastp_quality_metrics = {}
    
    try:
        with open(json_file_path, 'r', encoding='utf-8') as file:
            fastp_data = json.load(file)
            
        summary_data = fastp_data.get('summary', {})
        before_filtering_data = summary_data.get('before_filtering', {})
        
        # 提取并转换质量指标
        total_bases_gb = round(before_filtering_data.get('total_bases', 0) / 1e9, 2)
        q20_rate_percent = round(before_filtering_data.get('q20_rate', 0) * 100, 2)
        q30_rate_percent = round(before_filtering_data.get('q30_rate', 0) * 100, 2)
        gc_content_percent = round(before_filtering_data.get('gc_content', 0) * 100, 2)
        
        fastp_quality_metrics = {
            'total_bases': total_bases_gb,
            'q20_rate': q20_rate_percent,
            'q30_rate': q30_rate_percent, 
            'gc_content': gc_content_percent
        }
        
    except Exception as error:
        print(f"解析fastp JSON文件错误 {json_file_path}: {error}")
        # 设置默认值
        fastp_quality_metrics = {
            'total_bases': 0,
            'q20_rate': 0,
            'q30_rate': 0, 
            'gc_content': 0
        }
    
    return fastp_quality_metrics


def generate_json_qc_report(output_directory: str, sample_name: str, 
                           qc_variables: Dict[str, Any], 
                           fastp_quality_data: Dict[str, float], 
                           sample_id: str) -> None:
    """
    生成JSON格式的质控报告
    
    参数:
        output_directory: 输出目录
        sample_name: 样本名称
        qc_variables: 质控变量字典
        fastp_quality_data: fastp质量数据
        sample_id: 样本ID
    
    返回:
        None
    """
    
    # 构建JSON数据结构
    quality_control_report = {
        "sample_info": {
            "sample_name": sample_name,
            "sample_id": sample_id,
            "analysis_date": str(Path(__file__).parent.resolve())
        },
        "sequencing_quality": {
            "total_base_gb": fastp_quality_data['total_bases'],
            "q20_rate_percent": fastp_quality_data['q20_rate'],
            "q30_rate_percent": fastp_quality_data['q30_rate'],
            "gc_content_percent": fastp_quality_data['gc_content'],
            "total_clean_reads": qc_variables['Toatl_clean_read'],
            "mapping_rate_percent": qc_variables['Mapping_rate'],
            "duplication_rate_percent": qc_variables['Dup_rate']
        },
        "panel_region": {
            "target_capture_rate_percent": qc_variables['panelCapture_rate'],
            "average_depth": qc_variables['panelavedepth'],
            "coverage_1x_percent": qc_variables['panelCoverage'],
            "coverage_5x_percent": qc_variables['cov5_rate'],
            "coverage_10x_percent": qc_variables['panelcov10_rate'],
            "coverage_20x_percent": qc_variables['cov20_rate'],
            "coverage_30x_percent": qc_variables['cov30_rate'],
            "coverage_50x_percent": qc_variables['cov50_rate'],
            "coverage_100x_percent": qc_variables['panelcov100_rate'],
            "fload_80": round(qc_variables.get('fload80', 0), 2),
            "mapped_bases_gb": qc_variables['panelmap_base'],
            "mapped_reads_million": qc_variables['panelmap_read'],
            "target_bases_gb": qc_variables['paneltarget_base']
        }
    }
    
    # 写入JSON文件
    output_file_path = Path(output_directory) / f"{sample_name}.sample_information.json"
    with open(output_file_path, 'w', encoding='utf-8') as file:
        json.dump(quality_control_report, file, indent=4, ensure_ascii=False)
    
    print(f"JSON质控报告已生成: {output_file_path}")


def generate_txt_qc_report(output_directory: str, sample_name: str, 
                          qc_variables: Dict[str, Any], 
                          fastp_quality_data: Dict[str, float]) -> None:
    """
    生成TXT格式的质控报告
    
    参数:
        output_directory: 输出目录
        sample_name: 样本名称
        qc_variables: 质控变量字典
        fastp_quality_data: fastp质量数据
    
    返回:
        None
    """
    
    output_file_path = Path(output_directory) / f"{sample_name}.sample_information.txt"
    
    with open(output_file_path, 'w', encoding='utf-8') as file:
        # ==================== Panel区域信息 ====================
        file.write("panel information\n")
        file.write(f"Total base(G)\t{fastp_quality_data['total_bases']}\n")
        file.write(f"q20_rate(%)\t{fastp_quality_data['q20_rate']}\n")
        file.write(f"q30_rate(%)\t{fastp_quality_data['q30_rate']}\n")
        file.write(f"gc_content(%)\t{fastp_quality_data['gc_content']}\n")
        file.write(f"Total clean read\t{qc_variables['Toatl_clean_read']}\n")
        file.write(f"Mapping rate(%)\t{qc_variables['Mapping_rate']}\n")
        file.write(f"Duplication rate(%)\t{qc_variables['Dup_rate']}\n")
        file.write(f"Target capture rate(%)\t{qc_variables['panelCapture_rate']:.2f}\n")
        file.write(f"Average Depth on Target\t{qc_variables['panelavedepth']:.2f}\n")
        file.write(f"Coverage on Target(1X)\t{qc_variables['panelCoverage']:.2f}%\n")
        file.write(f"Coverage on Target(5X)\t{qc_variables['cov5_rate']:.4f}\n")
        file.write(f"Coverage on Target(10X)\t{qc_variables['panelcov10_rate']:.4f}\n")
        file.write(f"Coverage on Target(20X)\t{qc_variables['cov20_rate']:.4f}\n")
        file.write(f"Coverage on Target(30X)\t{qc_variables['cov30_rate']:.4f}\n")
        file.write(f"Coverage on Target(50X)\t{qc_variables['cov50_rate']:.4f}\n")
        file.write(f"Coverage on Target(100X)\t{qc_variables['panelcov100_rate']:.4f}\n")
        file.write(f"FLOAD 80 for Target\t{qc_variables.get('fload80', 0):.14f}\n")
    
    print(f"TXT质控报告已生成: {output_file_path}")


if __name__ == "__main__":
    # 命令行参数验证
    if len(sys.argv) != 7:
        print("用法: python qc_fload_analysis.py <sample> <bamfile> <dedup_metrics> <panelfile> <outdir> <sampleid>")
        print("参数说明:")
        print("  sample: 样本名称")
        print("  bamfile: BAM文件路径")
        print("  dedup_metrics: 重复标记质控结果文件路径")
        print("  panelfile: panel区域文件路径")
        print("  outdir: 输出目录")
        print("  sampleid: 样本ID")
        sys.exit(1)
    
    # 提取命令行参数
    sample_arg = sys.argv[1]
    bamfile_arg = sys.argv[2]
    dedup_metrics_arg = sys.argv[3]
    panelfile_arg = sys.argv[4]
    outdir_arg = sys.argv[5]     
    sampleid_arg = sys.argv[6]   
    
    # 执行质控分析
    qc_fload_analysis(sample_arg, bamfile_arg, dedup_metrics_arg, panelfile_arg, 
                     outdir_arg, sampleid_arg)