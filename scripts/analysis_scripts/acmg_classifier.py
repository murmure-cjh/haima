#!/usr/bin/env python3
"""
ACMG变异分类器模块
基于PathCall.pl逻辑实现ACMG分类规则
"""

import os
import re
from typing import Dict, Set, Any

# =============================================================================
# 默认配置 - 当没有配置文件或配置文件出错时使用
# =============================================================================

DEFAULT_CONFIG = {
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


class ACMGClassifier:
    """ACMG变异分类器，基于PathCall.pl逻辑实现"""
    
    def __init__(self, config_path: str = None, config_dict: Dict = None):
        """
        初始化ACMG分类器
        
        Args:
            config_path: YAML配置文件路径
            config_dict: 直接传入的配置字典
        """
        # 加载配置
        self.config = self._load_config(config_path, config_dict)
        
        # 获取ACMG分类器配置
        acmg_config = self.config.get('acmg_classifier', {})
        
        # 使用配置中的路径和参数
        self.db_paths = {
            "clinvar": acmg_config.get('clinvar', ''),
            "bic_brca1": acmg_config.get('bic_brca1', ''),
            "bic_brca2": acmg_config.get('bic_brca2', ''),
            "tp53": acmg_config.get('tp53', ''),
        }
        
        # 分类参数
        self.high_frequency_threshold = acmg_config.get('high_frequency_threshold', 0.05)
        self.splicing_distance_threshold = acmg_config.get('splicing_distance_threshold', 2)
        self.cds_end_threshold = acmg_config.get('cds_end_threshold', 50)
        
        # 氨基酸三字母到单字母代码映射
        self.codon_map = {
            "Ala": "A", "Cys": "C", "Asp": "D", "Glu": "E",
            "Phe": "F", "Gly": "G", "His": "H", "Ile": "I",
            "Lys": "K", "Leu": "L", "Met": "M", "Asn": "N",
            "Pro": "P", "Gln": "Q", "Arg": "R", "Ser": "S",
            "Thr": "T", "Val": "V", "Trp": "W", "Tyr": "Y",
            "Ter": "X"
        }
        
        # 验证路径
        self._validate_paths()
        
        # 加载ACMG数据库
        self.load_acmg_databases()

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
            return self._merge_configs(DEFAULT_CONFIG, config_dict)
        
        if config_path:
            try:
                import yaml
                if os.path.exists(config_path):
                    with open(config_path, 'r') as f:
                        loaded_config = yaml.safe_load(f)
                    return self._merge_configs(DEFAULT_CONFIG, loaded_config)
                else:
                    print(f"警告: ACMG配置文件不存在: {config_path}，使用默认配置")
            except Exception as e:
                print(f"警告: 加载ACMG配置文件失败: {e}，使用默认配置")
        
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
        required_paths = ["clinvar", "bic_brca1", "bic_brca2", "tp53"]
        
        missing_paths = []
        for path_key in required_paths:
            path = self.db_paths[path_key]
            if not path:
                missing_paths.append(f"{path_key}: 路径为空")
            elif not os.path.exists(path):
                missing_paths.append(f"{path_key}: {path}")
        
        if missing_paths:
            print("警告: 以下ACMG数据库文件不存在:")
            for missing in missing_paths:
                print(f"  - {missing}")

    def load_acmg_databases(self):
        """加载ACMG分类所需的所有数据库"""
        print("加载ACMG数据库...")
        self.acmg_dbs = {
            'clinvar': {},
            'bic': {},
            'tp53': {},
            'pathogenic_variants': set(),
            'benign_variants': set(),
            'pathogenic_aa_changes': {},
            'benign_aa_changes': {}
        }
        
        # 依次加载各个数据库
        if os.path.exists(self.db_paths["clinvar"]):
            print(f"加载ClinVar数据库: {self.db_paths['clinvar']}")
            self._load_clinvar_database()
        
        if os.path.exists(self.db_paths["bic_brca1"]):
            print(f"加载BIC BRCA1数据库: {self.db_paths['bic_brca1']}")
            self._load_bic_database(self.db_paths["bic_brca1"], "BRCA1")
        
        if os.path.exists(self.db_paths["bic_brca2"]):
            print(f"加载BIC BRCA2数据库: {self.db_paths['bic_brca2']}")
            self._load_bic_database(self.db_paths["bic_brca2"], "BRCA2")
        
        if os.path.exists(self.db_paths["tp53"]):
            print(f"加载TP53数据库: {self.db_paths['tp53']}")
            self._load_tp53_database()

    def _load_clinvar_database(self):
        """加载ClinVar数据库并解析致病/良性变异"""
        with open(self.db_paths["clinvar"], 'r') as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 5:
                    # 构建变异唯一标识
                    key = f"{parts[0]},{parts[1]},{parts[2]},{parts[3]}"
                    self.acmg_dbs['clinvar'][key] = parts[4]
                    
                    # 分类致病性和良性变异
                    if 'Pathogenic' in parts[4] or 'Likely pathogenic' in parts[4]:
                        self.acmg_dbs['pathogenic_variants'].add(key)
                        self._record_aa_change(parts, self.acmg_dbs['pathogenic_aa_changes'])
                    
                    if 'Benign' in parts[4] or 'Likely benign' in parts[4]:
                        self.acmg_dbs['benign_variants'].add(key)
                        self._record_aa_change(parts, self.acmg_dbs['benign_aa_changes'])

    def _load_bic_database(self, file_path: str, gene: str):
        """
        加载BIC数据库
        
        Args:
            file_path: 数据库文件路径
            gene: 基因名称(BRCA1或BRCA2)
        """
        with open(file_path, 'r') as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 5:
                    key = f"{parts[0]},{parts[1]},{parts[2]},{parts[3]}"
                    interpretation = f"{parts[3]}{parts[4]}"
                    self.acmg_dbs['bic'][key] = interpretation
                    
                    # 分类致病性和良性变异
                    if 'yes;' in parts[3] or 'Class 5' in parts[4]:
                        self.acmg_dbs['pathogenic_variants'].add(key)
                        if len(parts) >= 3 and parts[2].startswith('p.') and 'Ter' not in parts[2]:
                            aa_change = self._translate_aa_change(parts[2].replace(';', ''))
                            if gene not in self.acmg_dbs['pathogenic_aa_changes']:
                                self.acmg_dbs['pathogenic_aa_changes'][gene] = set()
                            self.acmg_dbs['pathogenic_aa_changes'][gene].add(aa_change)
                    
                    if 'no;' in parts[3] or 'Class 1' in parts[4]:
                        self.acmg_dbs['benign_variants'].add(key)

    def _load_tp53_database(self):
        """加载TP53数据库"""
        with open(self.db_paths["tp53"], 'r') as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 4:
                    key = f"{parts[0]},{parts[1]},{parts[2]},{parts[3]}"
                    self.acmg_dbs['tp53'][key] = parts[3]
                    
                    # 分类致病性和良性变异
                    if 'non-functional' in parts[3]:
                        self.acmg_dbs['pathogenic_variants'].add(key)
                        if len(parts) >= 3 and parts[2].startswith('p.') and 'Ter' not in parts[2]:
                            if 'TP53' not in self.acmg_dbs['pathogenic_aa_changes']:
                                self.acmg_dbs['pathogenic_aa_changes']['TP53'] = set()
                            self.acmg_dbs['pathogenic_aa_changes']['TP53'].add(parts[2])
                    
                    if ';functional;' in parts[3] or 'supertrans' in parts[3]:
                        self.acmg_dbs['benign_variants'].add(key)

    def _record_aa_change(self, parts: list, aa_dict: dict):
        """
        记录氨基酸变化到指定字典
        
        Args:
            parts: 数据库行分割后的列表
            aa_dict: 目标氨基酸变化字典
        """
        if len(parts) >= 5 and parts[2].startswith('p.') and 'Ter' not in parts[2] and parts[1]:
            aa_change = self._translate_aa_change(parts[2])
            gene = parts[1]
            if gene not in aa_dict:
                aa_dict[gene] = set()
            aa_dict[gene].add(aa_change)

    def _translate_aa_change(self, aa_change: str) -> str:
        """
        将氨基酸变化从三字母代码转换为单字母代码
        
        Args:
            aa_change: 三字母代码的氨基酸变化字符串
            
        Returns:
            单字母代码的氨基酸变化字符串
        """
        for three_letter, one_letter in self.codon_map.items():
            aa_change = aa_change.replace(three_letter, one_letter)
        return aa_change

    def classify_variant(self, variant: Dict, report_genes: Set) -> str:
        """
        对变异进行ACMG分类
        
        Args:
            variant: 变异信息字典
            report_genes: 报告基因集合
            
        Returns:
            ACMG分类结果字符串
        """
        # 构建变异唯一标识
        chrom = variant.get('Chr', '')
        pos = str(variant.get('Start', ''))
        ref = variant.get('ref', '')
        alt = variant.get('Alt', '')
        mutation_key = f"{chrom},{pos},{ref},{alt}"
        
        gene = variant.get('Gene.refGene', '')
        func = variant.get('Func.refGene', '')
        exonic_func = variant.get('ExonicFunc.refGene', '')
        
        # 检查是否在报告基因列表中
        in_report_genes = gene in report_genes
        
        # 检查是否在致病数据库中
        in_pathogenic_db = mutation_key in self.acmg_dbs['pathogenic_variants']
        in_benign_db = mutation_key in self.acmg_dbs['benign_variants']
        
        # 检查频率信息
        has_high_freq = self._has_high_frequency(variant)
        
        # ACMG分类逻辑（基于PathCall.pl）
        
        # 1. 高频率变异 -> 良性
        if has_high_freq:
            return "Benign"
        
        # 2. 不在报告基因中但存在于致病数据库 -> VGUS(Report pathogenic)
        if not in_report_genes and in_pathogenic_db:
            return "VGUS(Report pathogenic)"
        
        # 3. 不在报告基因中 -> VGUS
        if not in_report_genes:
            return "VGUS"
        
        # 4. 同义突变
        if exonic_func == "synonymous SNV":
            if in_pathogenic_db:
                return "Reported pathogenic"
            else:
                return "Likely Benign"
        
        # 5. 非编码区域变异
        if self._is_non_coding_region(func):
            if in_pathogenic_db:
                return "Reported pathogenic"
            else:
                return "Likely Benign"
        
        # 6. 剪接区域变异
        if func == "splicing":
            return self._classify_splicing_variant(variant, mutation_key, in_pathogenic_db, in_benign_db)
        
        # 7. 外显子区域变异
        if func == "exonic":
            return self._classify_exonic_variant(variant, mutation_key, gene, exonic_func, in_pathogenic_db, in_benign_db)
        
        # 默认分类
        return "VUS"

    def _has_high_frequency(self, variant: Dict) -> bool:
        """
        检查变异是否在人群中具有高频率
        
        Args:
            variant: 变异信息字典
            
        Returns:
            如果任何频率字段值 >= 阈值，返回True，否则False
        """
        freq_fields = [
            '1000g2015aug_eur', 'esp6500siv2_ea', 'ExAC_ALL', 'gnomAD_genome_ALL'
        ]
        
        for field in freq_fields:
            freq = variant.get(field, '0')
            if freq != '.' and freq != '0' and float(freq) >= self.high_frequency_threshold:
                return True
        return False

    def _is_non_coding_region(self, func: str) -> bool:
        """
        检查变异是否位于非编码区域
        
        Args:
            func: 变异功能区域
            
        Returns:
            如果是非编码区域返回True，否则False
        """
        non_coding_regions = ['downstream', 'upstream', 'ncRNA', 'intergenic', 'intronic', 'UTR5', 'UTR3']
        return func in non_coding_regions

    def _classify_splicing_variant(self, variant: Dict, mutation_key: str, in_pathogenic_db: bool, in_benign_db: bool) -> str:
        """
        分类剪接区域变异
        
        Args:
            variant: 变异信息字典
            mutation_key: 变异唯一标识
            in_pathogenic_db: 是否在致病数据库中
            in_benign_db: 是否在良性数据库中
            
        Returns:
            剪接区域变异的ACMG分类
        """
        splicing_distance = self._get_splicing_distance(variant)
        freq_1kg = variant.get('1000g2015aug_eur', '0')
        
        if splicing_distance <= self.splicing_distance_threshold:
            # 剪接位点变异（距离<=阈值）
            if in_pathogenic_db:
                # 检查ClinVar冲突
                clinvar_sig = self.acmg_dbs['clinvar'].get(mutation_key, '')
                if 'benign' in clinvar_sig.lower() or 'uncertain' in clinvar_sig.lower():
                    return "CLinVar Conflict report"
                elif freq_1kg != '.' and freq_1kg != '0':
                    return "1000genome exists"
                else:
                    return "Pathogenic"
            elif in_benign_db:
                return "Report Benign"
            elif freq_1kg != '.' and freq_1kg != '0':
                return "1000genome exists"
            else:
                return "Likely pathogenic"
        else:
            # 内含子变异（距离>阈值）
            if in_pathogenic_db:
                return "Reported pathogenic"
            elif in_benign_db:
                return "Likely Benign"
            else:
                return "VUS"

    def _classify_exonic_variant(self, variant: Dict, mutation_key: str, gene: str, 
                               exonic_func: str, in_pathogenic_db: bool, in_benign_db: bool) -> str:
        """
        分类外显子区域变异
        
        Args:
            variant: 变异信息字典
            mutation_key: 变异唯一标识
            gene: 基因名称
            exonic_func: 外显子功能类型
            in_pathogenic_db: 是否在致病数据库中
            in_benign_db: 是否在良性数据库中
            
        Returns:
            外显子区域变异的ACMG分类
        """
        # 截断突变（无义突变、移码突变）
        if exonic_func in ["stopgain", "stoploss", "frameshift insertion", "frameshift deletion"]:
            return self._classify_truncating_variant(variant, mutation_key, in_pathogenic_db, in_benign_db)
        
        # 错义突变和非移码突变
        cds_position = self._get_cds_position(variant)
        aa_change = variant.get('AAChange.refGene', '')
        
        # 检查起始密码子变异
        if cds_position <= 3:
            return self._classify_initiation_codon_variant(variant, mutation_key, in_pathogenic_db, in_benign_db)
        
        # 常规错义突变
        if in_pathogenic_db:
            return "Reported pathogenic"
        # 检查相同氨基酸变化
        elif self._has_same_aa_change(gene, aa_change):
            return "Same AA_change"
        # 检查软件预测有害
        elif self._is_probably_damaging(variant):
            return "VUS(Porobably Damaging)"
        elif in_benign_db:
            return "Likely Benign"
        else:
            return "VUS"

    def _classify_truncating_variant(self, variant: Dict, mutation_key: str, 
                                   in_pathogenic_db: bool, in_benign_db: bool) -> str:
        """
        分类截断突变
        
        Args:
            variant: 变异信息字典
            mutation_key: 变异唯一标识
            in_pathogenic_db: 是否在致病数据库中
            in_benign_db: 是否在良性数据库中
            
        Returns:
            截断变异的ACMG分类
        """
        cds_position = self._get_cds_position(variant)
        transcript = variant.get('Transcript', '')
        freq_1kg = variant.get('1000g2015aug_eur', '0')
        
        if in_pathogenic_db:
            # 检查ClinVar冲突
            clinvar_sig = self.acmg_dbs['clinvar'].get(mutation_key, '')
            if 'benign' in clinvar_sig.lower() or 'uncertain' in clinvar_sig.lower():
                return "ClinVar Conflict report"
            # 检查是否在CDS末端
            elif self._is_cds_end(transcript, cds_position):
                return "CDS end"
            elif freq_1kg != '.' and freq_1kg != '0':
                return "1000genome exists"
            else:
                return "Pathogenic"
        elif in_benign_db:
            return "Report Benign"
        else:
            # 新发截断突变
            if self._is_cds_end(transcript, cds_position):
                return "CDS end"
            elif freq_1kg != '.' and freq_1kg != '0':
                return "1000genome exists"
            else:
                return "Likely pathogenic"

    def _classify_initiation_codon_variant(self, variant: Dict, mutation_key: str, 
                                         in_pathogenic_db: bool, in_benign_db: bool) -> str:
        """
        分类起始密码子变异
        
        Args:
            variant: 变异信息字典
            mutation_key: 变异唯一标识
            in_pathogenic_db: 是否在致病数据库中
            in_benign_db: 是否在良性数据库中
            
        Returns:
            起始密码子变异的ACMG分类
        """
        freq_1kg = variant.get('1000g2015aug_eur', '0')
        
        if in_pathogenic_db:
            clinvar_sig = self.acmg_dbs['clinvar'].get(mutation_key, '')
            if 'benign' in clinvar_sig.lower() or 'uncertain' in clinvar_sig.lower():
                return "ClinVar Conflict report"
            elif freq_1kg != '.' and freq_1kg != '0':
                return "1000genome exists"
            else:
                return "Pathogenic"
        elif in_benign_db:
            return "Report Benign"
        elif freq_1kg != '.' and freq_1kg != '0':
            return "1000genome exists"
        else:
            return "Likely pathogenic"

    def _get_splicing_distance(self, variant: Dict) -> int:
        """
        获取剪接距离
        
        Args:
            variant: 变异信息字典
            
        Returns:
            剪接距离，默认返回10
        """
        info = variant.get('Otherinfo', '') or variant.get('Info', '')
        if 'splice' in info.lower():
            match = re.search(r'splice.*?[+-]?(\d+)', info)
            if match:
                return int(match.group(1))
        return 10  # 默认返回较大的距离

    def _get_cds_position(self, variant: Dict) -> int:
        """
        获取CDS位置
        
        Args:
            variant: 变异信息字典
            
        Returns:
            CDS位置，如果无法获取返回0
        """
        c_hgvs = variant.get('cHGVS', '')
        if c_hgvs:
            match = re.search(r'c\.(\d+)', c_hgvs)
            if match:
                return int(match.group(1))
        return 0

    def _is_cds_end(self, transcript: str, cds_position: int) -> bool:
        """
        检查是否在CDS末端
        
        Args:
            transcript: 转录本名称
            cds_position: CDS位置
            
        Returns:
            如果在CDS末端返回True，否则False
        """
        return cds_position > 0 and cds_position >= self.cds_end_threshold

    def _has_same_aa_change(self, gene: str, aa_change: str) -> bool:
        """
        检查是否有相同的氨基酸变化
        
        Args:
            gene: 基因名称
            aa_change: 氨基酸变化
            
        Returns:
            如果存在相同的氨基酸变化返回True，否则False
        """
        if gene in self.acmg_dbs['pathogenic_aa_changes']:
            translated_aa = self._translate_aa_change(aa_change)
            return translated_aa in self.acmg_dbs['pathogenic_aa_changes'][gene]
        return False

    def _is_probably_damaging(self, variant: Dict) -> bool:
        """
        检查软件预测是否为有害
        
        Args:
            variant: 变异信息字典
            
        Returns:
            如果SIFT和PolyPhen2都预测为有害返回True，否则False
        """
        sift_pred = variant.get('SIFT_pred', '')
        polyphen2_hdiv_pred = variant.get('Polyphen2_HDIV_pred', '')
        polyphen2_hvar_pred = variant.get('Polyphen2_HVAR_pred', '')
        
        return (sift_pred == 'D' and 
                polyphen2_hdiv_pred == 'D' and 
                polyphen2_hvar_pred == 'D')