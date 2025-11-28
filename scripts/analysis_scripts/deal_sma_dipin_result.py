#!/usr/bin/env python3
"""
合并后的DIPIN和SMA结果解析脚本
功能：解析地中海贫血和脊髓性肌肉萎缩症的检测结果，生成标准化的CSV报告
"""

import os
import sys
import re
import argparse
import csv
from pathlib import Path
from typing import List, Optional, Tuple

class GeneticResultParser:
    """遗传检测结果解析器"""
    
    def __init__(self, default_race: str = "Asian"):
        self.results = []
        self.default_race = default_race
    
    def parse_dipin_results(self, alpha_file: str, beta_file: str, other_file: str, 
                          sample: str, output_dir: str, race: str = None) -> str:
        """
        解析DIPIN（地中海贫血）检测结果
        
        Args:
            alpha_file: α地中海贫血文件路径
            beta_file: β地中海贫血文件路径
            other_file: 其他血红蛋白变异文件路径
            sample: 样本编号
            output_dir: 输出目录
            race: 种族信息
            
        Returns:
            生成的CSV文件路径
        """
        # 检查输入文件是否存在
        for file_path, file_name in [(alpha_file, "α地中海贫血"), 
                                   (beta_file, "β地中海贫血"), 
                                   (other_file, "其他血红蛋白变异")]:
            if not os.path.exists(file_path):
                print(f"错误: {file_name}文件不存在: {file_path}")
                sys.exit(1)
        
        # 创建输出目录 - 在输出目录下创建dipin子目录
        dipin_output_dir = os.path.join(output_dir, "dipin")
        os.makedirs(dipin_output_dir, exist_ok=True)
        
        output_file = os.path.join(dipin_output_dir, f"{sample}_dipin_result.csv")
        
        # 写入CSV文件头
        with open(output_file, 'w', newline='', encoding='utf-8') as f:
            writer = csv.writer(f)
            writer.writerow(["样本编号", "疾病类型", "突变基因", "结果", "染色体区间", "变异类型"])
        
        # 解析α地中海贫血文件
        self._parse_alpha_thalassemia(alpha_file, sample, output_file)
        
        # 解析β地中海贫血文件
        self._parse_beta_thalassemia(beta_file, sample, output_file)
        
        # 解析其他血红蛋白变异文件
        self._parse_other_hemoglobin(other_file, sample, output_file)
        
        # 根据种族信息修正结果
        if race:
            self._correct_for_ethnicity(output_file, race, sample)
        
        print(f"DIPIN结果生成成功: {output_file}")
        return output_file
    
    def _parse_alpha_thalassemia(self, file_path: str, sample: str, output_file: str):
        """解析α地中海贫血结果"""
        try:
            with open(file_path, 'r', encoding='utf-8') as f:
                for line in f:
                    line = line.strip()
                    if not line:
                        continue
                    
                    # 按制表符分割字段
                    fields = line.split('\t')
                    if len(fields) < 3:
                        continue
                    
                    chromosome_region = fields[0]
                    result = fields[1]
                    variant_type = fields[2]
                    
                    # 确定基因类型
                    if "HBA1" in line:
                        gene = "HBA1"
                    elif "HBA2" in line:
                        gene = "HBA2"
                    else:
                        gene = "-"
                    
                    # 写入结果
                    with open(output_file, 'a', newline='', encoding='utf-8') as out_f:
                        writer = csv.writer(out_f)
                        writer.writerow([sample, "α地中海贫血", gene, result, chromosome_region, variant_type])
                        
        except Exception as e:
            print(f"解析α地中海贫血文件错误: {e}")
            sys.exit(1)
    
    def _parse_beta_thalassemia(self, file_path: str, sample: str, output_file: str):
        """解析β地中海贫血结果"""
        beta_genes = ["HBE1", "HBE2", "HBG1", "HBG2", "HBD", "HBB"]
        
        try:
            with open(file_path, 'r', encoding='utf-8') as f:
                for line in f:
                    line = line.strip()
                    if not line:
                        continue
                    
                    fields = line.split('\t')
                    if len(fields) < 3:
                        continue
                    
                    chromosome_region = fields[0]
                    result = fields[1]
                    variant_type = fields[2]
                    
                    # 确定基因类型
                    gene = "-"
                    for beta_gene in beta_genes:
                        if beta_gene in line:
                            gene = beta_gene
                            break
                    
                    # 写入结果
                    with open(output_file, 'a', newline='', encoding='utf-8') as out_f:
                        writer = csv.writer(out_f)
                        writer.writerow([sample, "β地中海贫血", gene, result, chromosome_region, variant_type])
                        
        except Exception as e:
            print(f"解析β地中海贫血文件错误: {e}")
            sys.exit(1)
    
    def _parse_other_hemoglobin(self, file_path: str, sample: str, output_file: str):
        """解析其他血红蛋白变异结果"""
        try:
            with open(file_path, 'r', encoding='utf-8') as f:
                for line in f:
                    line = line.strip()
                    if not line:
                        continue
                    
                    fields = line.split('\t')
                    if len(fields) < 3:
                        continue
                    
                    chromosome_region = fields[0]
                    result = fields[1]
                    variant_type = fields[2]
                    
                    # 写入结果（基因列始终为空）
                    with open(output_file, 'a', newline='', encoding='utf-8') as out_f:
                        writer = csv.writer(out_f)
                        writer.writerow([sample, "其他血红蛋白变异", "-", result, chromosome_region, variant_type])
                        
        except Exception as e:
            print(f"解析其他血红蛋白变异文件错误: {e}")
            sys.exit(1)
    
    def _correct_for_ethnicity(self, result_file: str, race: str, sample: str):
        """根据种族信息修正结果"""
        try:
            # 检查结果文件中是否包含"Caucasian HPFH"
            with open(result_file, 'r', encoding='utf-8') as f:
                content = f.read()
            
            if "Caucasian HPFH" not in content:
                return
            
            # 检查种族信息是否为Asian
            if race == "Asian":
                # 替换Caucasian HPFH为SEA HPFH*
                with open(result_file, 'r', encoding='utf-8') as f:
                    content = f.read()
                
                content = content.replace("Caucasian HPFH", "SEA HPFH*")
                
                with open(result_file, 'w', encoding='utf-8') as f:
                    f.write(content)
                print(f"根据种族信息({race})修正了Caucasian HPFH为SEA HPFH*")
                    
        except Exception as e:
            print(f"种族信息修正错误: {e}")
    
    def parse_sma_results(self, sma_file: str, sample: str, output_dir: str) -> str:
        """
        解析SMA（脊髓性肌肉萎缩症）检测结果
        
        Args:
            sma_file: SMA结果文件路径
            sample: 样本编号
            output_dir: 输出目录
            
        Returns:
            生成的CSV文件路径
        """
        # 检查SMA文件是否存在
        if not os.path.exists(sma_file):
            print(f"错误: SMA结果文件不存在: {sma_file}")
            sys.exit(1)
        
        # 创建输出目录 - 在输出目录下创建sma子目录
        sma_output_dir = os.path.join(output_dir, "sma")
        os.makedirs(sma_output_dir, exist_ok=True)
        
        output_file = os.path.join(sma_output_dir, f"{sample}_SMA_result.csv")
        
        try:
            # 读取SMA文件内容
            with open(sma_file, 'r', encoding='utf-8') as f:
                content = f.read().strip()
            
            # 解析exon信息和基因型
            if '|' in content:
                parts = content.split('|')
                if len(parts) >= 2:
                    exon_info = parts[0].strip()
                    zo_info = parts[1].strip()
                    
                    # 处理exon信息，提取最后两个元素
                    exon_parts = exon_info.split('_')
                    if len(exon_parts) >= 2:
                        processed_exon = f"{exon_parts[-2]}_{exon_parts[-1]}"
                    else:
                        processed_exon = exon_info
                else:
                    processed_exon = "未知"
                    zo_info = "未知"
            else:
                processed_exon = "未知"
                zo_info = "未知"
            
            # 写入CSV文件
            with open(output_file, 'w', newline='', encoding='utf-8') as f:
                writer = csv.writer(f)
                # 写入文件头
                writer.writerow(["样本编号", "基因", "染色体位置", "转录本基因亚区", 
                               "MAF(人群频率)", "核苷酸/氨基酸", "基因型", "疾病(OMIM/遗传方式)"])
                # 写入数据行
                writer.writerow([sample, "SMN1", "/", "NM_000344.4 EX7", "/", 
                               processed_exon, zo_info, "脊髓性肌肉萎缩症(253550/AR)"])
            
            print(f"SMA结果生成成功: {output_file}")
            return output_file
            
        except Exception as e:
            print(f"解析SMA文件错误: {e}")
            sys.exit(1)

def main():
    """主函数"""
    parser = argparse.ArgumentParser(
        description="解析DIPIN和SMA检测结果，生成标准化报告",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
使用示例:
  python genetic_parser.py dipin --alpha alpha.txt --beta beta.txt --other other.txt --sample SAMPLE001 --output ./result --race Asian
  python genetic_parser.py sma --sma smacacnv.workout --sample SAMPLE002 --output ./result
  python genetic_parser.py all --alpha alpha.txt --beta beta.txt --other other.txt --sma smacacnv.workout --sample SAMPLE003 --output ./result --race Asian
        """
    )
    
    parser.add_argument("mode", choices=["dipin", "sma", "all"], 
                       help="运行模式: dipin-仅解析DIPIN, sma-仅解析SMA, all-解析全部")
    
    # DIPIN相关参数
    parser.add_argument("--alpha", help="α地中海贫血文件路径")
    parser.add_argument("--beta", help="β地中海贫血文件路径")
    parser.add_argument("--other", help="其他血红蛋白变异文件路径")
    
    # SMA相关参数
    parser.add_argument("--sma", help="SMA结果文件路径")
    
    # 通用参数
    parser.add_argument("--sample", required=True, help="样本编号")
    parser.add_argument("--output", required=True, help="输出目录路径")
    parser.add_argument("--race", default="Asian", help="种族信息(默认: Asian)")
    
    args = parser.parse_args()
    
    # 参数验证
    if args.mode in ["dipin", "all"]:
        if not args.alpha or not args.beta or not args.other:
            print("错误: DIPIN模式需要 --alpha, --beta 和 --other 参数")
            sys.exit(1)
    
    if args.mode in ["sma", "all"]:
        if not args.sma:
            print("错误: SMA模式需要 --sma 参数")
            sys.exit(1)
    
    parser = GeneticResultParser(default_race=args.race)
    
    try:
        if args.mode in ["dipin", "all"]:
            parser.parse_dipin_results(
                args.alpha, args.beta, args.other, 
                args.sample, args.output, args.race
            )
        
        if args.mode in ["sma", "all"]:
            parser.parse_sma_results(args.sma, args.sample, args.output)
            
        print("所有处理完成!")
        
    except Exception as e:
        print(f"处理过程中发生错误: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()