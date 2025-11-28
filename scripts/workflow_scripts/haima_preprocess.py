#!/usr/bin/env python3
"""
WES数据预处理脚本
功能：读取样本信息CSV，下载原始fq文件，生成snakemake配置文件
作者：chenjh
日期：2025.11

主要功能：
1. 命令行参数解析
2. 前置条件检查
3. 输入文件验证
4. 样本信息处理
5. 文件下载管理
6. 配置文件生成
7. 结果统计和报告
"""

import pandas as pd
import subprocess
import os
import yaml
import sys
import argparse
from pathlib import Path
from tqdm import tqdm
import time


class WESPreprocessor:
    """WES数据预处理主类，负责协调整个预处理流程"""
    
    def __init__(self):
        """初始化预处理器，设置统计信息和无效样本记录"""
        self.stats = {
            'total_samples': 0,
            'valid_samples': 0,
            'invalid_samples': 0,
            'successful_downloads': 0,
            'skipped_existing': 0,
            'failed_downloads': 0,
            'start_time': None,
            'end_time': None
        }
        self.invalid_samples_info = {}  # 存储无效样本及其原因
    
    def setup_argparse(self):
        """设置和解析命令行参数
        
        Returns:
            argparse.Namespace: 解析后的命令行参数
        """
        parser = argparse.ArgumentParser(
            description='WES数据预处理：下载fq文件并生成snakemake配置文件',
            formatter_class=argparse.RawDescriptionHelpFormatter,
            epilog=f'''
使用示例:
  python {os.path.basename(__file__)} -i input.csv -o config.yaml
  python {os.path.basename(__file__)} -i input.csv -o config.yaml -d /path/to/raw_data --dry-run
            '''
        )
        
        # 必需参数
        parser.add_argument('-i', '--input', required=True, 
                          help='输入CSV文件路径，包含样本信息')
        parser.add_argument('-o', '--output', required=True,
                          help='输出YAML配置文件路径，用于snakemake流程')
        
        # 可选参数
        parser.add_argument('-d', '--download-dir', 
                          default='/haplox/users/chenjh/haima/snakemake/raw_data/',
                          help='fq文件下载目录 (默认: /haplox/users/chenjh/haima/snakemake/raw_data/)')
        parser.add_argument('--cos-prefix', default='cos://sz-hapseq',
                          help='COS路径前缀 (默认: cos://sz-hapseq)')
        parser.add_argument('--dry-run', action='store_true',
                          help='模拟运行，不实际下载文件，用于测试')
        parser.add_argument('--force', action='store_true',
                          help='强制重新下载已存在的文件')
        parser.add_argument('--invalid-output', 
                          help='无效样本输出文件路径')
        parser.add_argument('--overwrite', action='store_true', default=True,
                          help='覆盖已存在的配置文件 (默认: True)')
        
        return parser.parse_args()
    
    
    def check_prerequisites(self):
        """检查运行所需的前置条件
        
        Returns:
            bool: 所有前置条件是否满足
        """
        print("🔍 检查前置条件...")
        
        # 检查coscli工具是否可用
        try:
            result = subprocess.run(['coscli', '--version'], 
                                  capture_output=True, text=True)
            if result.returncode == 0:
                print("✅ coscli工具可用")
                return True
            else:
                print("❌ coscli工具不可用")
                return False
        except Exception as e:
            print(f"❌ coscli工具检查失败: {e}")
            return False    


    def validate_input_file(self, csv_file):
        """验证输入CSV文件的完整性和格式
        
        Args:
            csv_file (str): 输入CSV文件路径
            
        Returns:
            bool: 文件是否有效
        """
        if not os.path.exists(csv_file):
            print(f"❌ 输入文件不存在: {csv_file}")
            return False
        
        try:
            df = pd.read_csv(csv_file)
            required_columns = ['Sample_ID', 'Lib_number', 'Raw_Path_R1', 'Raw_Path_R2', 'Information']
            missing_columns = [col for col in required_columns if col not in df.columns]
            
            if missing_columns:
                print(f"❌ CSV文件缺少必要列: {missing_columns}")
                return False
            
            # 检查Lib_number是否唯一
            if df['Lib_number'].duplicated().any():
                duplicates = df[df['Lib_number'].duplicated()]['Lib_number'].unique()
                print(f"⚠️  警告: 发现重复的Lib_number: {list(duplicates)}")
            
            self.stats['total_samples'] = len(df)
            print(f"✅ 输入文件验证通过，共 {len(df)} 个样本")
            return True
            
        except Exception as e:
            print(f"❌ CSV文件读取失败: {e}")
            return False
  
    
    def extract_gender_from_information(self, information):
        """从Information字段中提取性别信息
        
        Args:
            information (str): Information字段内容
            
        Returns:
            str: 提取的性别信息 ('Male', 'Female', 或 'Unknown')
        """
        if not isinstance(information, str):
            return "Unknown"
        
        if "男" in information:
            return "Male"
        elif "女" in information:
            return "Female"
        else:
            return "Unknown"
    
    def get_sample_location(self, lib_number):
        """根据Lib_number前缀确定样本采集地点
        
        Args:
            lib_number (str): 文库编号
            
        Returns:
            str: 样本地点 ('JX', 'SZ', 或 'Unknown')
        """
        if not isinstance(lib_number, str):
            return "Unknown"
        
        if lib_number.startswith('JX'):
            return 'JX'
        elif lib_number.startswith('SZ'):
            return 'SZ'
        else:
            return 'Unknown'
        
    
    def validate_sample(self, row):
        """验证单个样本数据的有效性
        
        Args:
            row (pd.Series): 样本数据行
            
        Returns:
            bool: 样本是否有效
        """
        lib_number = row['Lib_number']
        
        # 检查Lib_number类型
        if not isinstance(lib_number, str):
            self.invalid_samples_info[lib_number] = f'Sample name {lib_number} not valid string'
            return False
        
        # 检查样本地点
        sample_location = self.get_sample_location(lib_number)
        if sample_location == 'Unknown':
            self.invalid_samples_info[lib_number] = f'Sample location unknown from Lib_number: {lib_number}'
            return False
        
        return True
    
    def download_file(self, cos_path, local_path, dry_run=False, force=False):
        """下载单个文件从COS到本地
        
        Args:
            cos_path (str): COS源文件路径
            local_path (str): 本地目标路径
            dry_run (bool): 是否模拟运行
            force (bool): 是否强制重新下载
            
        Returns:
            str: 下载结果状态 ('success', 'failed', 'skipped', 'simulated')
        """
        if os.path.exists(local_path) and not force:
            return 'skipped'
        
        if dry_run:
            print(f"   📋 模拟下载: {cos_path} -> {local_path}")
            return 'simulated'
        
        try:
            os.makedirs(os.path.dirname(local_path), exist_ok=True)
            
            # 使用简单的coscli命令，不指定日志路径
            cmd = ['coscli', 'cp', cos_path, local_path]
            
            print(f"   📝 执行命令: {' '.join(cmd)}")
            
            result = subprocess.run(cmd, capture_output=True, text=True)
            
            if result.returncode == 0:
                return 'success'
            else:
                print(f"   ❌ 下载失败: {result.stderr}")
                return 'failed'
                
        except Exception as e:
            print(f"   ❌ 下载异常: {e}")
            return 'failed'
    
    def process_samples(self, args):
        """处理所有样本数据，包括验证、下载和配置生成
        
        Args:
            args (argparse.Namespace): 命令行参数
            
        Returns:
            dict: 生成的配置数据
        """
        print("\n📊 开始处理样本...")
        
        # 读取CSV文件
        df = pd.read_csv(args.input)
        config_data = {}
        
        # 提取项目Sed_ID
        if 'Sed_ID' in df.columns and len(df) > 0:
            sed_id = df.iloc[0]['Sed_ID']
            config_data['Sed_ID'] = sed_id
            print(f"📋 项目Sed_ID: {sed_id}")
        else:
            config_data['Sed_ID'] = 'Unknown'
            print("⚠️  警告: 未找到Sed_ID列或CSV为空，设置Sed_ID为Unknown")
        
        # 使用进度条处理每个样本
        valid_samples = 0
        with tqdm(total=len(df), desc="样本处理进度", unit="样本") as pbar:
            for index, row in df.iterrows():
                sample_id = row['Sample_ID']
                lib_number = row['Lib_number']
                
                # 更新进度条描述
                pbar.set_description(f"处理 {lib_number}")
                
                # 验证样本有效性
                if not self.validate_sample(row):
                    pbar.update(1)
                    continue
                
                # 提取样本元数据
                gender = self.extract_gender_from_information(row.get('Information', ''))
                sample_location = self.get_sample_location(lib_number)
                
                # 构建文件路径
                raw_path_r1 = row['Raw_Path_R1']
                raw_path_r2 = row['Raw_Path_R2']
                cos_r1 = f"{args.cos_prefix}{raw_path_r1}"
                cos_r2 = f"{args.cos_prefix}{raw_path_r2}"
                
                # 生成本地文件路径
                r1_filename = os.path.basename(raw_path_r1)
                r2_filename = os.path.basename(raw_path_r2)
                local_r1_path = os.path.join(args.download_dir, r1_filename)
                local_r2_path = os.path.join(args.download_dir, r2_filename)
                
                # 下载R1和R2文件
                r1_result = self.download_file(cos_r1, local_r1_path, 
                                             args.dry_run, args.force)
                r2_result = self.download_file(cos_r2, local_r2_path, 
                                             args.dry_run, args.force)
                
                # 更新下载统计
                if r1_result == 'success' and r2_result == 'success':
                    self.stats['successful_downloads'] += 1
                elif r1_result == 'skipped' and r2_result == 'skipped':
                    self.stats['skipped_existing'] += 1
                else:
                    self.stats['failed_downloads'] += 1
                
                # 构建样本配置数据
                config_data[lib_number] = {
                    'lib_number': lib_number,
                    'Panel': row.get('Panel', ''),
                    'Sed_ID': row['Sed_ID'],
                    'Sample_ID': sample_id,
                    'Name': row['Name'],
                    'Local_R1': local_r1_path,
                    'Local_R2': local_r2_path,
                    'raw_path_r1': raw_path_r1,
                    "raw_path_r2": raw_path_r2,
                    'haima_project': row.get('haima_project', ''),                    
                    'Information': row.get('Information', ''),
                    'Gender': gender,  # 添加性别信息
                    'Sample_location': sample_location,  # 添加样本地点信息
                    'Sample_source': row.get('Sample_source', ''),
                    'Cancer_species': row.get('Cancer_species', ''),
                    'DNA_type': row.get('DNA_type', ''),
                    'Lane': row.get('Lane', ''),
                    'Adapter': row.get('Adapter', ''),
                    'Index1': row.get('Index1', ''),
                    'Index2': row.get('Index2', ''),
                    'Estimated_yield': row.get('Estimated_yield', ''),
                    'Raw_Reads_Num(M)': row.get('Raw_Reads_Num(M)', ''),
                    'Raw_Yield(G)': row.get('Raw_Yield(G)', ''),
                    'Raw_Q20(%)': row.get('Raw_Q20(%)', ''),
                    'Raw_Q30(%)': row.get('Raw_Q30(%)', ''),
                    'Raw_GC(%)': row.get('Raw_GC(%)', ''),
                    'Clean_Reads_Num(M)': row.get('Clean_Reads_Num(M)', ''),
                    'Clean_Yield(G)': row.get('Clean_Yield(G)', ''),
                    'Clean_Q20(%)': row.get('Clean_Q20(%)', ''),
                    'Clean_Q30(%)': row.get('Clean_Q30(%)', ''),
                    'Clean_GC(%)': row.get('Clean_GC(%)', ''),
                    'Effective(%)': row.get('Effective(%)', ''),
                    'Duplication_Rate(%)': row.get('Duplication_Rate(%)', ''),
                    'Data_ID': row.get('Data_ID', ''),
                    'Race': row.get('Race', ''),
                    'Download_Status': 'success' if (r1_result in ['success', 'skipped'] and 
                                                   r2_result in ['success', 'skipped']) else 'failed',
                    'Validation_Status': 'valid'
                }
                
                valid_samples += 1
                pbar.update(1)
        
        # 更新统计信息
        self.stats['valid_samples'] = valid_samples
        self.stats['invalid_samples'] = len(self.invalid_samples_info)
        
        return config_data
    
    def save_invalid_samples(self, df, output_file):
        """保存无效样本信息到CSV文件
        
        Args:
            df (pd.DataFrame): 原始数据框
            output_file (str): 输出文件路径
        """
        if not output_file or not self.invalid_samples_info:
            return
        
        try:
            # 筛选无效样本
            invalid_samples = []
            for index, row in df.iterrows():
                lib_number = row['Lib_number']
                if lib_number in self.invalid_samples_info:
                    invalid_row = row.copy()
                    invalid_row['invalid_reason'] = self.invalid_samples_info[lib_number]
                    invalid_samples.append(invalid_row)
            
            # 保存无效样本信息
            if invalid_samples:
                invalid_df = pd.DataFrame(invalid_samples)
                os.makedirs(os.path.dirname(output_file), exist_ok=True)
                invalid_df.to_csv(output_file, index=False)
                print(f"✅ 无效样本信息已保存: {output_file}")
                
        except Exception as e:
            print(f"❌ 保存无效样本信息失败: {e}")
    
    def save_config(self, config_data, output_file, overwrite=True):
        """保存配置数据到YAML文件
        
        Args:
            config_data (dict): 配置数据字典
            output_file (str): 输出文件路径
            overwrite (bool): 是否覆盖已存在文件
            
        Returns:
            bool: 保存是否成功
        """
        try:
            # 确保输出目录存在
            os.makedirs(os.path.dirname(output_file), exist_ok=True)
            
            # 检查文件是否存在
            if os.path.exists(output_file):
                if overwrite:
                    print(f"⚠️  配置文件已存在，将覆盖: {output_file}")
                else:
                    print(f"❌ 配置文件已存在且不允许覆盖: {output_file}")
                    return False
            
            # 写入YAML文件
            with open(output_file, 'w', encoding='utf-8') as f:
                yaml.dump(config_data, f, default_flow_style=False, 
                         allow_unicode=True, indent=2, sort_keys=False)
            
            print(f"✅ 配置文件已保存: {output_file}")
            return True
            
        except Exception as e:
            print(f"❌ 配置文件保存失败: {e}")
            return False
    
    def print_summary(self):
        """打印处理结果摘要"""
        print("\n" + "="*50)
        print("📈 处理摘要")
        print("="*50)
        print(f"总样本数: {self.stats['total_samples']}")
        print(f"有效样本: {self.stats['valid_samples']}")
        print(f"无效样本: {self.stats['invalid_samples']}")
        print(f"成功下载: {self.stats['successful_downloads']}")
        print(f"跳过已存在: {self.stats['skipped_existing']}")
        print(f"下载失败: {self.stats['failed_downloads']}")
        
        # 打印无效样本详情
        if self.invalid_samples_info:
            print("\n📋 无效样本详情:")
            for sample, reason in self.invalid_samples_info.items():
                print(f"  {sample}: {reason}")
        
        # 计算并显示处理时间
        if self.stats['start_time'] and self.stats['end_time']:
            duration = self.stats['end_time'] - self.stats['start_time']
            print(f"处理时间: {duration:.2f} 秒")
        
        print("="*50)
    
    def run(self):
        """主运行函数，协调整个预处理流程"""
        self.stats['start_time'] = time.time()
        
        print("🚀 WES数据预处理脚本启动")
        print("="*50)
        
        # 解析命令行参数
        args = self.setup_argparse()
        
        # 显示运行参数
        print("📋 运行参数:")
        print(f"  输入文件: {args.input}")
        print(f"  输出配置: {args.output}")
        print(f"  下载目录: {args.download_dir}")
        print(f"  COS前缀: {args.cos_prefix}")
        print(f"  模拟运行: {args.dry_run}")
        print(f"  强制下载: {args.force}")
        print(f"  覆盖配置: {args.overwrite}")
        if args.invalid_output:
            print(f"  无效样本输出: {args.invalid_output}")
        
        # 检查前置条件
        if not self.check_prerequisites():
            sys.exit(1)
        
        # 验证输入文件
        if not self.validate_input_file(args.input):
            sys.exit(1)
        
        # 读取原始数据用于无效样本记录
        df = pd.read_csv(args.input)
        
        # 处理所有样本
        config_data = self.process_samples(args)
        
        # 保存配置文件
        if not self.save_config(config_data, args.output, args.overwrite):
            sys.exit(1)
        
        # 保存无效样本信息
        invalid_output = args.invalid_output or f"invalid_samples_{os.path.basename(args.input)}"
        self.save_invalid_samples(df, invalid_output)
        
        self.stats['end_time'] = time.time()
        
        # 打印处理摘要
        self.print_summary()
        
        print("✅ 处理完成!")


def main():
    """主函数入口"""
    preprocessor = WESPreprocessor()
    preprocessor.run()


if __name__ == "__main__":
    main()