import os
import json
import glob
import pandas as pd
from pathlib import Path
import logging

def get_sed_id(batch_id):
    """根据 batch_id 获取 sed_id (前四段)"""
    parts = batch_id.split('_')
    if len(parts) >= 4:
        return "_".join(parts[:4])
    return batch_id

def generate_final_report(batch_id, base_work_dir):
    """
    功能：合并 QC JSON 结果到原始 CSV，生成最终报告。
    参数：
      batch_id: 批次号
      base_work_dir: snakemake 工作根目录
    返回：
      成功则返回生成的 CSV 绝对路径，失败返回 None
    """
    try:
        work_path = Path(base_work_dir)
        sed_id = get_sed_id(batch_id)
        
        logging.info(f"[Merge Module] Processing batch: {batch_id}, sed_id: {sed_id}")

        # 1. 定义路径
        json_search_path = work_path / "results" / sed_id
        raw_csv_path = work_path / "sample_info/raw_haima_csv" / f"{batch_id}.csv"
        
        # 2. 检查原始CSV是否存在
        if not raw_csv_path.exists():
            logging.error(f"[Merge Module] Raw CSV not found: {raw_csv_path}")
            return None

        # 3. 查找并解析 JSON 文件
        # 递归查找所有 .sample_information.json
        json_files = glob.glob(str(json_search_path / "**" / "*.sample_information.json"), recursive=True)
        
        if not json_files:
            logging.warning(f"[Merge Module] No JSON files found in {json_search_path}")
            # 如果找不到 JSON，根据需求决定是返回 None 还是仅返回原始 CSV
            # 这里选择返回 None，表示合并这一步“没做成”
            return None
            
        json_data_list = []
        for json_file in json_files:
            try:
                with open(json_file, 'r', encoding='utf-8') as f:
                    data = json.load(f)
                    # 提取关联键
                    sample_name = data.get("sample_info", {}).get("sample_name")
                    
                    if sample_name:
                        # 构造一行数据，Key 用 Lib_number 以便 merge
                        row = {"Lib_number": sample_name}
                        # 扁平化合并 sequencing_quality 和 depth_coverage
                        row.update(data.get("sequencing_quality", {}))
                        row.update(data.get("depth_coverage", {}))
                        row.update(data.get("panel_region", {}))
                        json_data_list.append(row)
            except Exception as e:
                logging.error(f"[Merge Module] Error reading {json_file}: {e}")

        if not json_data_list:
            logging.warning("[Merge Module] valid data extracted from JSONs is empty.")
            return None

        # 转为 DataFrame
        json_df = pd.DataFrame(json_data_list)
        
        # 4. 读取原始 CSV
        input_csv_df = pd.read_csv(raw_csv_path)
        
        if "Lib_number" not in input_csv_df.columns:
            logging.error(f"[Merge Module] 'Lib_number' column missing in {raw_csv_path}")
            return None

        # 5. 执行合并 (Left Join, 以原始 CSV 为主)
        final_df = pd.merge(input_csv_df, json_df, on="Lib_number", how="left")

        # 6. 保存结果
        # 输出文件名建议加上后缀，避免直接覆盖原始文件导致不可逆错误
        output_filename = f"{batch_id}_qc_report.csv"
        output_path = work_path / "results" / sed_id / output_filename
        
        # 保存为 UTF-8-SIG (防止 Excel 打开乱码)，保留3位小数
        final_df.to_csv(output_path, index=False, encoding="utf-8-sig", float_format='%.3f')
        
        logging.info(f"[Merge Module] Report generated: {output_path}")
        return str(output_path)

    except Exception as e:
        logging.error(f"[Merge Module] Critical error: {e}")
        return None