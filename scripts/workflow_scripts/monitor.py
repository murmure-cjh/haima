# a script to monitor cos://sz-hapseq/rawfq/JX_health/Nova/sample_info_merge folder for new csv file
import os
import subprocess
import datetime
import time
import logging
import sys

# 添加send_mail模块路径
sys.path.append('/haplox/users/chenjh/haima/snakemake/')
import send_mail
import merge_qc
# add log
logging.basicConfig(filename='/haplox/users/chenjh/haima/snakemake/logs/monitor.log', level=logging.INFO, 
                    format='%(asctime)s - %(levelname)s - %(message)s')

# track csv files - 只监控sample_info_merge文件夹
monitor_folder = 'cos://sz-hapseq/rawfq/JX_health/Nova/sample_info_merge/'
local_folder = '/haplox/users/chenjh/haima/snakemake/sample_info/raw_haima_csv/'
snakemake_work_dir = '/haplox/users/chenjh/haima/snakemake/'

# 任务队列
task_queue = []
currently_running = None
completed_tasks = set()  # 记录已完成的任务

# 确保本地目录存在
os.makedirs(local_folder, exist_ok=True)
os.makedirs('/haplox/users/chenjh/haima/snakemake/sample_info/snakemake_sample_yaml/', exist_ok=True)

# 获取已处理的文件列表
processed_files = set()
for file in os.listdir(local_folder):
    if '.csv' in file:
        processed_files.add(file)

def is_snakemake_running():
    """更精确的snakemake运行状态检查"""
    try:
        # 方法1：检查特定工作目录的snakemake进程
        result = subprocess.run([
            'pgrep', '-f', f'snakemake.*{snakemake_work_dir}'
        ], capture_output=True, text=True)
        
        # 方法2：检查锁文件（更可靠）
        lock_dir = os.path.join(snakemake_work_dir, '.snakemake/locking')
        if os.path.exists(lock_dir) and os.listdir(lock_dir):
            return True
            
        # 方法3：检查进程返回值
        if result.returncode == 0 and result.stdout.strip():
            return True
            
        return False
    except Exception as e:
        logging.error(f"Error checking snakemake status: {e}")
        return False  # 改为返回False，避免过度保守

def download_file(filename, folder=local_folder):
    try:
        logging.info(f"Starting download of file: {filename}")
        # 清理文件名
        filename = filename.split('|')[0].replace(' ', '')
        # 构建完整路径
        full_path = 'cos://sz-hapseq/' + filename
        # 获取真实文件名
        file_name_real = full_path.split('/')[-1]
        # 下载文件到本地文件夹
        os.system('coscli cp ' + full_path + ' ' + folder + '/' + file_name_real)
        logging.info(f"Downloaded file: {filename}")
        return file_name_real
    except Exception as e:
        logging.error(f"Error downloading file: {filename} : {e}")
        return None

def send_start_notification(batch_id, csv_file_path):
    """发送任务开始通知，附带CSV文件作为附件"""
    try:
        subject = send_mail.start_email_subject.format(batch_id)
        content = send_mail.start_email_content.format(batch_id)
        # 检查CSV文件是否存在
        if os.path.exists(csv_file_path):
            send_mail.send_email(send_mail.maintainer_receviers, subject, content, csv_file_path)
            logging.info(f"Sent start notification for batch: {batch_id} with attachment: {csv_file_path}")
        else:
            logging.warning(f"CSV file not found for attachment: {csv_file_path}, sending email without attachment")
            send_mail.send_email(send_mail.maintainer_receviers, subject, content)
    except Exception as e:
        logging.error(f"Error sending start notification for batch {batch_id}: {e}")

def send_completion_notification(batch_id, attachment_path=None):
    """
    发送任务完成通知
    :param batch_id: 批次号
    :param attachment_path: (可选) 附件的绝对路径，通常是 Final Report CSV
    """
    try:
        subject = send_mail.end_email_subject.format(batch_id)
        content = send_mail.end_email_content.format(batch_id)
        
        # 判断是否存在附件路径且文件有效
        if attachment_path and os.path.exists(attachment_path):
            logging.info(f"Found attachment for completion email: {attachment_path}")
            send_mail.send_email(send_mail.maintainer_receviers, subject, content, attachment_path)
            logging.info(f"Sent completion notification for batch: {batch_id} WITH attachment.")
        else:
            # 如果没有生成报告，或者路径不存在，则发送不带附件的邮件
            if attachment_path:
                logging.warning(f"Attachment file not found: {attachment_path}")
            
            send_mail.send_email(send_mail.maintainer_receviers, subject, content)
            logging.info(f"Sent completion notification for batch: {batch_id} WITHOUT attachment.")
            
    except Exception as e:
        logging.error(f"Error sending completion notification for batch {batch_id}: {e}")

def run_snakemake_workflow(batch_id):
    try:
        logging.info(f"Starting snakemake workflow for batch: {batch_id}")
        
        # 构建CSV文件路径
        csv_file_path = os.path.join(local_folder, f"{batch_id}.csv")
        
        # 发送开始通知（附带CSV附件）
        send_start_notification(batch_id, csv_file_path)
        
        # 切换到工作目录
        original_dir = os.getcwd()
        os.chdir(snakemake_work_dir)
                
        # 执行数据预处理脚本
        preprocess_cmd = f'python /haplox/users/chenjh/haima/snakemake/scripts/workflow_scripts/haima_preprocess.py -i sample_info/raw_haima_csv/{batch_id}.csv -o sample_info/snakemake_sample_yaml/{batch_id}.yaml'
        preprocess_res = subprocess.run(preprocess_cmd, shell=True, executable='/bin/bash')
        if preprocess_res.returncode != 0:
            logging.error("Preprocess failed")
            return None       
        # 运行snakemake流程
        snakemake_cmd = f'snakemake --cores 36 -p --config sample_config="sample_info/snakemake_sample_yaml/{batch_id}.yaml" --rerun-incomplete >> logs/snakemake.log 2>&1'
        
        # 执行所有命令
        full_cmd = f"source /haplox/users/chenjh/miniforge3/bin/activate /haplox/users/chenjh/miniforge3/envs/snakemake && {snakemake_cmd}"
        proc = subprocess.Popen(full_cmd, shell=True, executable='/bin/bash')
        return proc
        
    except Exception as e:
        logging.error(f"Failed to start workflow: {e}")
        return None
# 假设你在 monitor.py 全局定义了 current_process 用于追踪运行中的 Popen 对象
# global current_process 

def process_queue():
    global currently_running, current_process, completed_tasks

    # 1. 检查当前是否有任务在运行
    if currently_running and current_process:
        # 轮询进程状态
        ret_code = current_process.poll()
        
        if ret_code is not None:  # 进程已结束
            if ret_code == 0:
                logging.info(f"Batch {currently_running} Snakemake workflow finished successfully.")
                
                try:
                    # 1. 调用 merge_qc 生成最终报告
                    logging.info(f"Starting QC merge for batch: {currently_running}")
                    # 传入 batch_id 和 snakemake 的工作根目录
                    final_report_path = merge_qc.generate_final_report(currently_running, snakemake_work_dir)
                    
                    if final_report_path:
                        logging.info(f"QC Merge successful. Report generated at: {final_report_path}")
                    else:
                        logging.warning(f"QC Merge returned None for batch: {currently_running}")

                    # 2. 发送完成邮件（带上刚刚生成的报告作为附件）
                    send_completion_notification(currently_running, attachment_path=final_report_path)
                    
                except Exception as e:
                    logging.error(f"Error during Merge/Email phase for batch {currently_running}: {e}")
                    # 即使合并出错，可能也需要发送一个不带附件的通知，或者报警
                    send_completion_notification(currently_running)
                # ================= 核心修改结束 =================
                
                completed_tasks.add(currently_running)
                
            else:
                # 任务失败逻辑
                logging.error(f"Batch {currently_running} failed with return code {ret_code}. Check snakemake logs.")
                # 可选：发送失败报警邮件
            
            # 重置状态，准备处理下一个任务
            currently_running = None
            current_process = None

    # 2. 如果当前空闲且队列有任务，启动新任务
    if currently_running is None and task_queue:
        next_batch = task_queue.pop(0)
        # 注意：run_snakemake_workflow 需要修改为返回 Popen 对象并赋值给 global current_process
        logging.info(f"Starting processing for batch: {next_batch}")
        currently_running = next_batch
        current_process = run_snakemake_workflow(next_batch)

# 监控新文件上传
def process_new_files(monitor_folder, download_folder=local_folder): 
    try:
        logging.info(f"Starting processing of new files in folder: {monitor_folder}")
        (status1, output1) = subprocess.getstatusoutput("coscli ls " + monitor_folder)
        # 提取所有csv文件
        csv_files = [line for line in output1.split('\n') if '.csv' in line]
        new_batches = []
        # 遍历每个文件
        for file in csv_files:
            clean_file = file.split('|')[0].split('/')[-1].replace(' ', '')
            # 如果文件未被处理
            if clean_file not in processed_files:
                downloaded_file = download_file(file, download_folder)
                if downloaded_file:
                    processed_files.add(clean_file)
                    # 获取批次ID（去掉.csv扩展名）
                    batch_id = clean_file.split('.csv')[0]
                    new_batches.append(batch_id)
        return new_batches
    except Exception as e:
        logging.error(f"Error processing new files in folder: {monitor_folder} : {e}")
        return []

if __name__ == "__main__":
    logging.info("Script started")
    while True:
        logging.info("Starting new iteration")
        
        # 检查是否有新文件
        new_batches = process_new_files(monitor_folder)
        
        # 如果有新文件，加入队列
        for batch_id in new_batches:
            if batch_id not in task_queue and batch_id != currently_running and batch_id not in completed_tasks:
                task_queue.append(batch_id)
                logging.info(f"Added batch to queue: {batch_id}. Queue length: {len(task_queue)}")
        
        # 处理队列
        process_queue()
        
        # 记录队列状态
        if task_queue:
            logging.info(f"Current queue: {task_queue}. Currently running: {currently_running}")
        else:
            logging.info("Queue is empty")
        
        # 每30分钟检查一次
        now = datetime.datetime.now()
        formatted_date = now.strftime("%Y-%m-%d %H:%M:%S")
        logging.info(f"Sleeping for 10 minutes at {formatted_date}")
        time.sleep(600)