#!/usr/bin/env python3
"""
海马WES数据处理脚本
集成ACMG分类器的变异注释和过滤流程
"""

import json
import os
import sys
import pandas as pd
import re
import yaml
from typing import Dict, List, Any, Set

# 导入ACMG分类器
from acmg_classifier import ACMGClassifier

# =============================================================================
# 默认配置 - 当没有配置文件或配置文件出错时使用
# =============================================================================

DEFAULT_CONFIG = {
    # 海马处理器数据库配置
    "haima_processor": {
        "database_root": "/haplox/users/chenjh/database/haima_file/",
        "hgmd": "/haplox/users/chenjh/database/haima_file/hgmd_pro_2018.3_hg19.vcf",
        "omim": "/haplox/users/chenjh/database/haima_file/morbidmap.txt",
        "pseudogene": "/haplox/users/chenjh/database/haima_file/hgnc_pseudogene.txt",
        "vaf_baseline": "/haplox/users/chenjh/database/haima_file/10936sample.snv.baseline",
        "gene_correction": "/haplox/users/chenjh/database/haima_file/gene_anno_db.tsv",
        "intron_anno": "/haplox/users/chenjh/database/haima_file/intron_anno_db.tsv",
        "report_genes": "/haplox/users/chenjh/database/haima_file/LRGClinVarHGNC_germline.txt",
    },
    
    # ACMG分类器配置
    "acmg_classifier": {
        "database_root": "/haplox/users/chenjh/database/haima_file/ACMG",
        "clinvar": "/haplox/users/chenjh/database/haima_file/ACMG/ClinVar_format.txt",
        "bic_brca1": "/haplox/users/chenjh/database/haima_file/ACMG/BICBRCA1_format.txt",
        "bic_brca2": "/haplox/users/chenjh/database/haima_file/ACMG/BICBRCA2_format.txt",
        "tp53": "/haplox/users/chenjh/database/haima_file/ACMG/TP53_format.txt",
        "high_frequency_threshold": 0.05,
        "splicing_distance_threshold": 2,
        "cds_end_threshold": 50,
    }
}


class HaimaProcessor:
    """海马WES数据处理器，集成变异注释、过滤和ACMG分类功能"""
    
    def __init__(self, config_path: str = None, config_dict: Dict = None):
        """
        初始化处理器
        
        Args:
            config_path: YAML配置文件路径
            config_dict: 直接传入的配置字典
        """
        # 加载配置
        self.config = self._load_config(config_path, config_dict)
        
        # 获取海马处理器配置
        haima_config = self.config.get('haima_processor', {})
        
        # 使用配置中的路径
        self.paths = {
            # HGMD/OMIM相关数据库
            "hgmd": haima_config.get('hgmd', ''),
            "omim": haima_config.get('omim', ''),
            "pseudogene": haima_config.get('pseudogene', ''),
            "vaf_baseline": haima_config.get('vaf_baseline', ''),
            
            # 矫正数据库
            "gene_correction": haima_config.get('gene_correction', ''),
            "intron_anno": haima_config.get('intron_anno', ''),
            
            # 报告基因列表
            "report_genes": haima_config.get('report_genes', ''),
        }
        
        # 验证必要路径是否存在
        self._validate_paths()
        
        # 加载所有数据库
        self.load_all_databases()
        
        # 初始化ACMG分类器
        self.acmg_classifier = ACMGClassifier(config_dict=self.config)

    def _load_config(self, config_path: str = None, config_dict: Dict = None) -> Dict:
        """
        加载配置，支持多种配置来源
        
        Args:
            config_path: YAML配置文件路径
            config_dict: 直接传入的配置字典
            
        Returns:
            配置字典
        """
        # 优先级: config_dict > config_path > 默认配置
        if config_dict:
            print("使用传入的配置字典")
            return self._merge_configs(DEFAULT_CONFIG, config_dict)
        
        if config_path:
            try:
                if os.path.exists(config_path):
                    print(f"加载配置文件: {config_path}")
                    with open(config_path, 'r') as f:
                        loaded_config = yaml.safe_load(f)
                    return self._merge_configs(DEFAULT_CONFIG, loaded_config)
                else:
                    print(f"警告: 配置文件不存在: {config_path}，使用默认配置")
            except Exception as e:
                print(f"警告: 加载配置文件失败: {e}，使用默认配置")
        
        print("使用默认配置")
        return DEFAULT_CONFIG.copy()

    def _merge_configs(self, default: Dict, custom: Dict) -> Dict:
        """
        深度合并配置字典
        
        Args:
            default: 默认配置
            custom: 自定义配置
            
        Returns:
            合并后的配置字典
        """
        result = default.copy()
        
        for key, value in custom.items():
            if (key in result and 
                isinstance(result[key], dict) and 
                isinstance(value, dict)):
                # 递归合并字典
                result[key] = self._merge_configs(result[key], value)
            else:
                # 直接覆盖
                result[key] = value
        
        return result

    def _validate_paths(self):
        """验证必要路径是否存在"""
        required_paths = [
            "hgmd", "omim", "pseudogene", "vaf_baseline",
            "gene_correction", "intron_anno", "report_genes"
        ]
        
        missing_paths = []
        for path_key in required_paths:
            path = self.paths[path_key]
            if not path:
                missing_paths.append(f"{path_key}: 路径为空")
            elif not os.path.exists(path):
                missing_paths.append(f"{path_key}: {path}")
        
        if missing_paths:
            print("警告: 以下数据库文件不存在:")
            for missing in missing_paths:
                print(f"  - {missing}")

    def load_all_databases(self):
        """加载所有必要的数据库"""
        print("加载数据库...")
        self.hgmd_db = self.load_hgmd_database()
        self.omim_db = self.load_omim_database()
        self.pseudogene_db = self.load_pseudogene_database()
        self.vaf_baseline_db = self.load_vaf_baseline_database()
        self.gene_correction_db = self.load_gene_correction_database()
        self.intron_anno_db = self.load_intron_anno_database()
        self.report_genes_db = self.load_report_genes_database()

    def load_hgmd_database(self) -> Dict:
        """
        加载HGMD数据库
        
        Returns:
            HGMD数据库字典，键为变异标识，值为分类
        """
        hgmd_db = {}
        if os.path.exists(self.paths["hgmd"]):
            print(f"加载HGMD数据库: {self.paths['hgmd']}")
            with open(self.paths["hgmd"], 'r') as f:
                for line in f:
                    if line.startswith('#'):
                        continue
                    parts = line.strip().split()
                    if len(parts) < 8:
                        continue
                    # 解析INFO字段获取致病性分类
                    info_parts = parts[7].split(';')
                    cl_part = info_parts[0].split('=')
                    if len(cl_part) == 2 and cl_part[1] in ["DM", "DM?"]:
                        key = f"{parts[0]}\t{parts[1]}\t{parts[3]}\t{parts[4]}"
                        hgmd_db[key] = cl_part[1]
        else:
            print(f"警告: HGMD数据库文件不存在: {self.paths['hgmd']}")
        return hgmd_db

    def load_omim_database(self) -> Dict:
        """
        加载OMIM数据库
        
        Returns:
            OMIM数据库字典，键为基因名，值为OMIM编号
        """
        omim_db = {}
        if os.path.exists(self.paths["omim"]):
            print(f"加载OMIM数据库: {self.paths['omim']}")
            with open(self.paths["omim"], 'r') as f:
                for line in f:
                    if line.startswith('#'):
                        continue
                    parts = line.strip().split('\t')
                    if len(parts) < 2:
                        continue
                    # 解析基因列表
                    genes = parts[1].split(', ')
                    for gene in genes:
                        if gene not in omim_db:
                            omim_db[gene] = parts[0]
                        else:
                            omim_db[gene] += f";{parts[0]}"
        else:
            print(f"警告: OMIM数据库文件不存在: {self.paths['omim']}")
        return omim_db

    def load_pseudogene_database(self) -> Set:
        """
        加载伪基因数据库
        
        Returns:
            伪基因名称集合
        """
        pseudogene_db = set()
        if os.path.exists(self.paths["pseudogene"]):
            print(f"加载伪基因数据库: {self.paths['pseudogene']}")
            with open(self.paths["pseudogene"], 'r') as f:
                next(f)  # 跳过标题行
                for line in f:
                    parts = line.strip().split('\t')
                    if len(parts) > 1:
                        pseudogene_db.add(parts[1])
        else:
            print(f"警告: 伪基因数据库文件不存在: {self.paths['pseudogene']}")
        return pseudogene_db

    def load_vaf_baseline_database(self) -> Dict:
        """
        加载VAF基线数据库
        
        Returns:
            VAF基线数据库字典，键为变异标识，值为频率
        """
        vaf_baseline_db = {}
        if os.path.exists(self.paths["vaf_baseline"]):
            print(f"加载VAF基线数据库: {self.paths['vaf_baseline']}")
            with open(self.paths["vaf_baseline"], 'r') as f:
                for line in f:
                    parts = line.strip().split('\t')
                    if len(parts) >= 9:
                        key = f"{parts[4]}\t{parts[5]}\t{parts[6]}\t{parts[7]}"
                        vaf_baseline_db[key] = parts[8]
        else:
            print(f"警告: VAF基线数据库文件不存在: {self.paths['vaf_baseline']}")
        return vaf_baseline_db

    def load_gene_correction_database(self) -> Dict:
        """
        加载基因矫正数据库
        
        Returns:
            基因矫正字典，键为旧基因名，值为新基因名
        """
        gene_correction_db = {}
        if os.path.exists(self.paths["gene_correction"]):
            print(f"加载基因矫正数据库: {self.paths['gene_correction']}")
            try:
                df = pd.read_csv(self.paths["gene_correction"], sep='\t')
                gene_correction_db = dict(zip(df['Gene_old'], df['Gene']))
            except Exception as e:
                print(f"加载基因矫正数据库错误: {e}")
        else:
            print(f"警告: 基因矫正数据库文件不存在: {self.paths['gene_correction']}")
        return gene_correction_db

    def load_intron_anno_database(self) -> Dict:
        """
        加载内含子注释数据库
        
        Returns:
            内含子注释数据库字典
        """
        intron_anno_db = {}
        if os.path.exists(self.paths["intron_anno"]):
            print(f"加载内含子注释数据库: {self.paths['intron_anno']}")
            try:
                df = pd.read_csv(self.paths["intron_anno"], sep='\t')
                if 'chr_pos_ref_mut' in df.columns:
                    for col in df.columns:
                        if col != 'chr_pos_ref_mut':
                            intron_anno_db[col] = dict(zip(df['chr_pos_ref_mut'], df[col]))
            except Exception as e:
                print(f"加载内含子注释数据库错误: {e}")
        else:
            print(f"警告: 内含子注释数据库文件不存在: {self.paths['intron_anno']}")
        return intron_anno_db

    def load_report_genes_database(self) -> Set:
        """
        加载报告基因列表
        
        Returns:
            报告基因名称集合
        """
        report_genes = set()
        if os.path.exists(self.paths["report_genes"]):
            print(f"加载报告基因列表: {self.paths['report_genes']}")
            with open(self.paths["report_genes"], 'r') as f:
                next(f)  # 跳过标题行
                for line in f:
                    parts = line.strip().split('\t')
                    if len(parts) > 0:
                        report_genes.add(parts[0])
        else:
            print(f"警告: 报告基因列表文件不存在: {self.paths['report_genes']}")
        return report_genes

    def calculate_depth_info(self, variant: Dict) -> Dict:
        """
        计算测序深度信息
        
        Args:
            variant: 变异信息字典
            
        Returns:
            深度信息字典，包含参考深度、变异深度、总深度和变异频率
        """
        depth_info = {
            'refDepth': 0,
            'altDepth': 0,
            'allDepth': 0,
            'Readsratio': 0
        }
        
        # 从AD字段解析深度信息
        if 'AD' in variant and variant['AD']:
            try:
                ad_parts = variant['AD'].strip().split('/')
                if len(ad_parts) >= 2:
                    ref_depth = int(ad_parts[0])
                    alt_depth = int(ad_parts[1])
                    total_depth = ref_depth + alt_depth
                    
                    depth_info.update({
                        'refDepth': ref_depth,
                        'altDepth': alt_depth,
                        'allDepth': total_depth,
                        'Readsratio': alt_depth / total_depth if total_depth > 0 else 0
                    })
            except (ValueError, IndexError):
                pass
        # 如果AD字段不可用，使用DP字段
        elif 'DP' in variant and variant['DP']:
            try:
                total_depth = int(variant['DP'])
                depth_info.update({
                    'allDepth': total_depth,
                    'refDepth': total_depth // 2,
                    'altDepth': total_depth // 2,
                    'Readsratio': 0.5
                })
            except ValueError:
                pass
        
        return depth_info

    def determine_zygosity_from_gt(self, gt_value: str) -> str:
        """
        根据GT字段值判断合子类型
        
        Args:
            gt_value: GT字段值，格式如 "0/1", "1/1", "0|1", "1|1" 等
            
        Returns:
            中文合子类型: '纯合', '杂合', 或 '半合'
        """
        if not gt_value or gt_value in ['.', '-', '']:
            return '纯合'  # 默认值
        
        # 标准化GT值，移除分隔符前后的空格
        gt_clean = gt_value.strip().replace(' ', '')
        
        # 常见的GT格式及其含义:
        # 0/0, 0|0: 纯合参考 (通常不记录为变异)
        # 0/1, 0|1, 1/0, 1|0: 杂合
        # 1/1, 1|1: 纯合变异
        # 0: 半合参考 (男性X/Y染色体)
        # 1: 半合变异 (男性X/Y染色体)
        # 0/2, 1/2: 复合杂合
        
        # 检查是否为纯合变异 (1/1, 1|1)
        if re.match(r'^1[/|]1$', gt_clean):
            return '纯合'
        
        # 检查是否为杂合 (0/1, 0|1, 1/0, 1|0)
        elif re.match(r'^[01][/|][01]$', gt_clean) and gt_clean not in ['0/0', '0|0']:
            return '杂合'
        
        # 检查是否为复合杂合 (涉及两个不同等位基因)
        elif re.match(r'^[0-9]+[/|][0-9]+$', gt_clean) and len(set(gt_clean.split('/'))) > 1:
            return '杂合'
        
        # 检查是否为半合 (单个等位基因)
        elif gt_clean in ['0', '1']:
            return '半合'
        
        # 默认情况
        else:
            return '纯合'

    def correct_gene_info(self, variant: Dict) -> Dict:
        """
        矫正基因名和转录本信息
        
        Args:
            variant: 原始变异信息字典
            
        Returns:
            矫正后的变异信息字典
        """
        corrected_variant = variant.copy()
        
        # 构建位置键用于查找
        chrom = variant.get('Chr', '')
        pos = str(variant.get('Start', ''))
        ref = variant.get('Vcf_mut', '').split(':')[-1].split('/')[0] if 'Vcf_mut' in variant else variant.get('ref', '')
        alt = variant.get('Alt', '')
        
        position_key = f"{chrom};{pos};{ref};{alt}"
        
        # 应用基因名矫正
        current_gene = variant.get('Gene.refGene', '')
        if current_gene in self.gene_correction_db:
            corrected_variant['Gene.refGene'] = self.gene_correction_db[current_gene]
        
        # 应用内含子注释矫正
        for field, mapping in self.intron_anno_db.items():
            if position_key in mapping:
                corrected_variant[field] = mapping[position_key]
        
        return corrected_variant

    def translate_mutation_type(self, mutation_type: str) -> str:
        """
        翻译突变类型为中文
        
        Args:
            mutation_type: 英文突变类型
            
        Returns:
            中文突变类型
        """
        translation_map = {
            'nonsynonymous SNV': '错义突变',
            'synonymous SNV': '同义突变',
            'frameshift deletion': '移码缺失',
            'frameshift insertion': '移码插入',
            'stopgain': '无义突变',
            'stoploss': '无义突变',
            'nonframeshift insertion': '非移码插入',
            'nonframeshift deletion': '非移码缺失',
            'intronic': '内含子突变',
            'splicing': '剪接突变',
            'exonic;splicing': '剪接突变'
        }
        return translation_map.get(mutation_type, mutation_type)

    def translate_classification(self, classification: str) -> str:
        """
        翻译ACMG分类为中文
        
        Args:
            classification: 英文ACMG分类
            
        Returns:
            中文ACMG分类
        """
        translation_map = {
            'Pathogenic': '致病突变',
            'Likely pathogenic': '疑似致病',
            'VUS': '意义不明',
            'VGUS': '意义不明',
            'CLinVar Conflict report': '证据冲突',
            '1000genome exists': '存在千人频率',
            'Reported pathogenic': '报道致病',
            'CDS end': 'CDS 截断',
            'Same AA_change': '相同氨基酸改变',
            'Porobably Damaging': '软件预测有害突变',
            'VGUS(Report pathogenic)': '意义不明(报道致病)',
            'VUS(Report pathogenic)': '意义不明(报道致病)',
            'Benign': '良性',
            'Likely Benign': '疑似良性',
            'VUS(Porobably Damaging)': '意义不明(软件预测有害)',
            'VGUS(Porobably Damaging)': '意义不明(软件预测有害)'
        }
        return translation_map.get(classification, classification)

    def translate_zygosity(self, zygosity: str) -> str:
        """
        翻译合子性为中文 - 现在基于GT字段判断
        
        Args:
            zygosity: GT字段值或合子性字符串
            
        Returns:
            中文合子性
        """
        # 如果传入的是GT字段值，使用新的判断逻辑
        if zygosity and ( '/' in zygosity or '|' in zygosity or zygosity in ['0', '1'] ):
            return self.determine_zygosity_from_gt(zygosity)
        
        # 原有的字符串映射作为后备
        if zygosity == 'Heterozygous':
            return '杂合'
        else:
            return '纯合'

    def filter_variants(self, variant: Dict, depth_info: Dict) -> bool:
        """
        过滤变异
        
        Args:
            variant: 变异信息字典
            depth_info: 深度信息字典
            
        Returns:
            如果变异通过过滤返回True，否则False
        """
        # 深度过滤 - 深度>=10
        if depth_info['allDepth'] < 10:
            return False
        
        # 转录本过滤
        transcript = variant.get('HGNC_ID', '') or variant.get('Transcript', '')
        if transcript == 'UNKNOWN':
            return False
        
        # 伪基因过滤
        gene = variant.get('Gene.refGene', '')
        if gene in self.pseudogene_db:
            return False
        
        # 功能区域过滤 - 跳过非编码区域
        func = variant.get('Func.refGene', '')
        if func in ['intergenic', 'ncRNA']:  #'downstream', 'upstream', 
            return False
        
        return True

    def process_variant(self, variant: Dict) -> Dict:
        """
        处理单个变异
        
        Args:
            variant: 原始变异信息字典
            
        Returns:
            处理后的变异信息字典
        """
        # 步骤1: 生成深度信息
        depth_info = self.calculate_depth_info(variant)
        
        # 步骤2: 矫正基因信息
        corrected_variant = self.correct_gene_info(variant)
        
        # 步骤3: 添加深度信息
        corrected_variant.update(depth_info)
        
        # 步骤4: 应用HGMD注释
        chrom = corrected_variant.get('Chr', '')
        pos = str(corrected_variant.get('Start', ''))
        ref = corrected_variant.get('Vcf_mut', '').split(':')[-1].split('/')[0] if 'Vcf_mut' in corrected_variant else corrected_variant.get('ref', '')
        alt = corrected_variant.get('Alt', '')
        
        db_key = f"{chrom.replace('chr', '')}\t{pos}\t{ref}\t{alt}" 
        corrected_variant['HGMD'] = self.hgmd_db.get(db_key, ".")
        
        # 步骤5: 应用OMIM注释
        gene = corrected_variant.get('Gene.refGene', '')
        corrected_variant['OMIM'] = self.omim_db.get(gene, ".")
        
        # 步骤6: 应用VAF基线数据
        corrected_variant['HapWES10000_AF'] = self.vaf_baseline_db.get(db_key, "0")
        
        # 步骤7: ACMG生物信息学分类
        acmg_classification = self.acmg_classifier.classify_variant(corrected_variant, self.report_genes_db)
        corrected_variant['ACMG_Classification'] = acmg_classification
        
        # 步骤8: 应用格式转换
        mutation_type = corrected_variant.get('ExonicFunc.refGene', '') or corrected_variant.get('Func.refGene', '')
        if mutation_type:
            corrected_variant['mutation_type_cn'] = self.translate_mutation_type(mutation_type)
        
        if acmg_classification:
            corrected_variant['classification_cn'] = self.translate_classification(acmg_classification)
        
        # 修改：优先使用GT字段判断合子性
        gt_value = corrected_variant.get('GT', '')
        if gt_value:
            corrected_variant['zygosity_cn'] = self.translate_zygosity(gt_value)
        else:
            # 如果没有GT字段，回退到原有的Zygosity字段
            zygosity = corrected_variant.get('Zygosity', '')
            corrected_variant['zygosity_cn'] = self.translate_zygosity(zygosity)
        
        return corrected_variant

    def process_json_data(self, input_file: str, output_file: str):
        """
        处理JSON格式的变异数据
        
        Args:
            input_file: 输入JSON文件路径
            output_file: 输出JSON文件路径
        """
        print(f"读取输入文件: {input_file}")
        with open(input_file, 'r') as f:
            data = json.load(f)
        
        # 处理变异数据
        print("处理变异数据...")
        processed_data = []
        filtered_count = 0
        
        for variant in data['data']:
            # 处理单个变异
            processed_variant = self.process_variant(variant)
            
            # 计算深度信息用于过滤
            depth_info = self.calculate_depth_info(variant)
            
            # 过滤变异
            if self.filter_variants(processed_variant, depth_info):
                processed_data.append(processed_variant)
            else:
                filtered_count += 1
        
        # 更新输出数据结构
        output_data = data.copy()
        output_data['data'] = processed_data
        
        # 添加新列到columns列表
        new_columns = [
            'HGMD', 'OMIM', 'HapWES10000_AF', 'ACMG_Classification',
            'refDepth', 'altDepth', 'allDepth', 'Readsratio',
            'mutation_type_cn', 'classification_cn', 'zygosity_cn'
        ]
        
        # 确保不重复添加列
        existing_columns = set(output_data['columns'])
        for col in new_columns:
            if col not in existing_columns:
                output_data['columns'].append(col)
        
        # 写入输出文件
        print(f"写入输出文件: {output_file}")
        with open(output_file, 'w') as f:
            json.dump(output_data, f, indent=2, ensure_ascii=False)
        
        print(f"处理完成!")
        print(f"总变异数: {len(data['data'])}")
        print(f"处理后变异数: {len(processed_data)}")
        print(f"过滤掉的变异: {filtered_count}")


def main():
    """主函数，处理命令行参数并启动处理流程"""
    if len(sys.argv) < 3:
        print("用法: python haimaresult.py -i input.json -o output.json [-c config.yaml]")
        print("示例: python haimaresult.py -i input.json -o output.json -c config.yaml")
        sys.exit(1)
    
    # 解析命令行参数
    input_file = None
    output_file = None
    config_file = None
    
    i = 1
    while i < len(sys.argv):
        if sys.argv[i] == "-i" and i + 1 < len(sys.argv):
            input_file = sys.argv[i + 1]
            i += 2
        elif sys.argv[i] == "-o" and i + 1 < len(sys.argv):
            output_file = sys.argv[i + 1]
            i += 2
        elif sys.argv[i] == "-c" and i + 1 < len(sys.argv):
            config_file = sys.argv[i + 1]
            i += 2
        else:
            i += 1
    
    if not all([input_file, output_file]):
        print("错误: 必须提供 -i 和 -o 参数")
        sys.exit(1)
    
    # 创建处理器并处理数据
    try:
        processor = HaimaProcessor(config_path=config_file)
        processor.process_json_data(input_file, output_file)
    except Exception as e:
        print(f"处理过程中发生错误: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()