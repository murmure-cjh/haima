# miniwes.smk - hg38/SMA/Dipin 分析模块
# 功能：执行hg38参考基因组的比对分析、SMA（脊髓性肌萎缩症）分析和Dipin分析
# 依赖：需要config.yaml配置文件和common_def.smk通用函数

# ==================== 配置和包含文件 ====================
configfile: "config.yaml"

# 包含通用函数定义
include: "common_def.smk"

# ==================== HG38比对分析模块 ====================

rule bwa_mem_hg38:
    """使用Sentieon优化的BWA MEM在hg38参考基因组上进行比对"""
    input:
        fq1 = rules.fastp.output.clean_fq1,
        fq2 = rules.fastp.output.clean_fq2,
        dirs_ready = rules.create_directories.output.directory_created
    output:
        bam = "results/{sed_id}/{sample}/hg38/{sample}.hg38.sort.bam"
    params:
        sentieon_dir = config['sentieon_dir'],
        sentieon_license = config['sentieon_license'],
        ref_hg38 = config['hg38']  
    log:
        "results/{sed_id}/{sample}/logs/bwa_mem_hg38.log"
    threads: config['threads']
    conda: "../env_yaml/xiaohaima_x03.yaml"
    shell:
        """
        echo "开始对样本 {wildcards.sample} 进行hg38 BWA比对..."
        
        export SENTIEON_LICENSE={params.sentieon_license}
        
        # 执行BWA比对和排序
        {params.sentieon_dir} bwa mem \
            -R "@RG\\tID:{wildcards.sample}\\tLB:{wildcards.sample}\\tPL:ILLUMINA\\tSM:{wildcards.sample}" \
            -t {threads} -k 32 -M {params.ref_hg38} \
            {input.fq1} {input.fq2} | \
        {params.sentieon_dir} util sort \
            -o {output.bam} -t {threads} --sam2bam -i -
        
        echo "{wildcards.sample}hg38 BWA比对完成，输出文件: {output.bam}"
        """

rule dedup_hg38:
    """在hg38参考基因组上进行PCR重复标记和去除"""
    input:
        bam = rules.bwa_mem_hg38.output.bam,
        dirs_ready = rules.create_directories.output.directory_created
    output:
        dedup_bam = "results/{sed_id}/{sample}/hg38/{sample}.hg38.dedup.bam",
        score = "results/{sed_id}/{sample}/tmp/{sample}.hg38.score.txt",
        metrics = "results/{sed_id}/{sample}/hg38/{sample}.hg38.dedup_metrics.txt"
    params:
        sentieon_dir = config['sentieon_dir'],
        sentieon_license = config['sentieon_license'],
        ref_hg38 = config['hg38']
    log:
        "results/{sed_id}/{sample}/logs/dedup_hg38.log"
    threads: config['threads']
    conda: "../env_yaml/xiaohaima_x03.yaml"
    shell:
        """
        echo "开始对样本 {wildcards.sample} 进行hg38去重处理..."
        
        export SENTIEON_LICENSE={params.sentieon_license}
        
        # 步骤1: 计算位点收集器得分信息
        {params.sentieon_dir} driver \
            -r {params.ref_hg38} -t {threads} -i {input.bam} \
            --algo LocusCollector --fun score_info {output.score}
        
        # 步骤2: 执行去重算法
        {params.sentieon_dir} driver \
            -t {threads} -i {input.bam} \
            --algo Dedup --rmdup --score_info {output.score} \
            --metrics {output.metrics} {output.dedup_bam} \
            2>{output.dedup_bam}.dup.log
        
        echo "{wildcards.sample}hg38去重处理完成，输出文件: {output.dedup_bam}"
        """

rule realign_bqsr_hg38:
    """在hg38参考基因组上执行局部重比对和碱基质量值重校准"""
    input:
        dedup_bam = rules.dedup_hg38.output.dedup_bam,
        dirs_ready = rules.create_directories.output.directory_created
    output:
        realigned_bam = "results/{sed_id}/{sample}/hg38/{sample}.hg38.realigned.bam",
        bqsr_bam = "results/{sed_id}/{sample}/hg38/{sample}_dipin.bqsr.bam",
        recal_table = "results/{sed_id}/{sample}/hg38/{sample}_dipin.recal_data.table"
    params:
        sentieon_dir = config['sentieon_dir'],
        sentieon_license = config['sentieon_license'],
        ref_hg38 = config['hg38'],
        dipin_bed_slop60 = config['dipin_bed_slop60'], 
        hg38_snps = config['hg38_snps'],  
        hg38_indels = config['hg38_indels'], 
        hg38_dbsnp = config['hg38_dbsnp']   
    log:
        realign = "results/{sed_id}/{sample}/logs/{sample}_hg38_realign.log",
        bqsr = "results/{sed_id}/{sample}/logs/{sample}_hg38_bqsr.log"
    threads: config['threads']
    conda: "../env_yaml/xiaohaima_x03.yaml"
    shell:
        """
        echo "开始对样本 {wildcards.sample} 进行hg38重比对和BQSR..."
        
        export SENTIEON_LICENSE={params.sentieon_license}
        
        # 步骤1: 局部重比对
        {params.sentieon_dir} driver \
            -r {params.ref_hg38} -t {threads} -i {input.dedup_bam} \
            --algo Realigner \
            -k {params.hg38_snps} -k {params.hg38_indels} -k {params.hg38_dbsnp} \
            --interval_list {params.dipin_bed_slop60} \
            {output.realigned_bam} \
            2> {log.realign}
        
        # 步骤2: 碱基质量值重校准
        {params.sentieon_dir} driver \
            -t {threads} -r {params.ref_hg38} -i {output.realigned_bam} \
            --interval {params.dipin_bed_slop60} \
            --algo QualCal \
            -k {params.hg38_snps} -k {params.hg38_indels} -k {params.hg38_dbsnp} \
            {output.recal_table}
        
        # 步骤3: 应用BQSR校正
        {params.sentieon_dir} driver \
            -t {threads} -r {params.ref_hg38} -i {output.realigned_bam} \
            --read_filter QualCalFilter,table={output.recal_table},indel=false \
            --interval {params.dipin_bed_slop60} \
            --algo ReadWriter \
            {output.bqsr_bam}
        
        echo "{wildcards.sample} hg38重比对和BQSR完成，输出文件: {output.bqsr_bam}"
        """

rule haplotyper_hg38:
    """在hg38参考基因组上进行变异检测"""
    input:
        bqsr_bam = rules.realign_bqsr_hg38.output.bqsr_bam,
        dirs_ready = rules.create_directories.output.directory_created
    output:
        hc_vcf = "results/{sed_id}/{sample}/hg38/{sample}.hc.vcf"
    params:
        sentieon_dir = config['sentieon_dir'],
        sentieon_license = config['sentieon_license'],
        ref_hg38 = config['hg38'],
        dipin_bed_slop60 = config['dipin_bed_slop60']
    log:
        "results/{sed_id}/{sample}/logs/hg38_Haplotyper.log"
    threads: config['threads']
    conda: "../env_yaml/xiaohaima_x03.yaml"
    shell:
        """
        echo "开始对样本 {wildcards.sample} 进行hg38变异检测..."
        
        export SENTIEON_LICENSE={params.sentieon_license}
        
        # 执行Haplotyper变异检测
        {params.sentieon_dir} driver \
            -r {params.ref_hg38} -t {threads} -i {input.bqsr_bam} \
            --interval {params.dipin_bed_slop60} \
            --algo Haplotyper \
            {output.hc_vcf} \
            2>> {log}
        
        echo "{wildcards.sample} hg38变异检测完成，输出文件: {output.hc_vcf}"
        """

# ==================== SMA分析模块 ====================

rule sma_analysis:
    """进行脊髓性肌萎缩症（SMA）相关分析"""
    input:
        dedup_bam = rules.dedup.output.dedup_bam,
        gender_file = rules.sample_gender_analysis.output.gender_file,
        dirs_ready = rules.create_directories.output.directory_created
    output:
        sma_csv = "results/{sed_id}/{sample}/sma/{sample}.smacacnv.csv",
        sma_workout = "results/{sed_id}/{sample}/sma/{sample}.smacacnv.workout.csv"
    params:
        dmd_bam_files = get_dmd_bam_files,  # 获取DMD相关BAM文件列表
        workdir = os.path.abspath("."),
        uid = os.getuid(),
        gid = os.getgid(),
        docker_image = config['docker']['image'],
        docker_volumes = config['docker']['volumes'],
        script_path = config['Docker_scripts']['deal_sma_workout']
    log:
        "results/{sed_id}/{sample}/logs/sma_analysis.log"
    conda: "../env_yaml/xiaohaima_x03.yaml"
    shell:
        """
        echo "开始对样本 {wildcards.sample} 进行SMA分析..."
        
        # 步骤1: 运行smaca进行SMA分析
        docker run --rm \
            --user {params.uid}:{params.gid} \
            -v {params.workdir}:/workspace \
            {params.docker_volumes} \
            {params.docker_image} \
            smaca --reference hg19 \
                  --output {output.sma_csv} \
                  {input.dedup_bam} {params.dmd_bam_files}
        
        # 步骤2: 运行SMA结果后处理脚本
        docker run --rm \
            --user {params.uid}:{params.gid} \
            -v {params.workdir}:/workspace \
            {params.docker_volumes} \
            {params.docker_image} \
            python /scripts/deal_sma_workout.py \
                {output.sma_csv} \
                {output.sma_workout} \
                {wildcards.sample}.dedup.bam
        
        echo "样本{wildcards.sample}————SMA分析完成，输出文件: {output.sma_csv}, {output.sma_workout}"
        """

# ==================== DIPIN分析模块 ====================
# DIPIN分析模块是基于一个单独的镜像，所以部分写法有所差异
# 相关的镜像是单独购买的nad的软件，软件原理详见说明文档
rule dipin_analysis:
    """执行Dipin分析（基于hg38的CNV分析）"""
    input:
        bqsr_bam = rules.realign_bqsr_hg38.output.bqsr_bam,
        hc_vcf = rules.haplotyper_hg38.output.hc_vcf,
        dirs_ready = rules.create_directories.output.directory_created
    output:
        dipin_alpha = "results/{sed_id}/{sample}/dipin/result/CNV/α地中海贫血.txt",
        dipin_beta = "results/{sed_id}/{sample}/dipin/result/CNV/β地中海贫血.txt",
        dipin_other = "results/{sed_id}/{sample}/dipin/result/CNV/其他血红蛋白变异.txt"
    params:
        docker_container = config['dipin_docker_container'],
        ref_hg38 = config['hg38'],
        hg38_annovar = config['hg38_annovar'], 
        hg38_snps = config['hg38_snps'],
        hg38_indels = config['hg38_indels'], 
        hg38_dbsnp = config['hg38_dbsnp'],
        workdir = lambda wildcards: os.path.abspath("."),
        uid = os.getuid(),
        gid = os.getgid(),
        docker_volumes = config['docker']['volumes']
    log:
        "results/{sed_id}/{sample}/logs/dipin_analysis.log"
    conda: "../env_yaml/xiaohaima_x03.yaml"
    shell:
        """
        echo "开始对样本 {wildcards.sample} 进行Dipin分析..."
        
        # 获取当前工作目录的绝对路径
        WORKDIR={params.workdir}

        # 修改输出目录权限
        docker run --rm \
            -v {params.workdir}:/workspace \
            {params.docker_volumes} \
            {params.docker_container} \
            chmod 777 -R /workspace/results/{wildcards.sed_id}/{wildcards.sample}/dipin/        
        
        # 步骤1: 运行Dipin分析docker容器
        docker run --rm \
            -v {params.workdir}:/workspace \
            {params.docker_volumes} \
            {params.docker_container} hgbp \
                -t {wildcards.sample} \
                -b /workspace/{input.bqsr_bam} \
                -v /workspace/{input.hc_vcf} \
                -w /workspace/results/{wildcards.sed_id}/{wildcards.sample}/dipin/ \
                -a {params.hg38_annovar} \
                -g {params.ref_hg38} \
                -s {params.hg38_snps} \
                -m {params.hg38_indels} \
                -d {params.hg38_dbsnp}
        
        # 修改输出目录权限
        docker run --rm \
            -v {params.workdir}:/workspace \
            {params.docker_volumes} \
            {params.docker_container} \
            chmod 777 -R /workspace/results/{wildcards.sed_id}/{wildcards.sample}/dipin/
        
        echo "样本{wildcards.sample}————Dipin分析完成"
        """

# ==================== SMA和DIPIN结果处理模块 ====================

rule deal_sma_dipin_result:
    """处理SMA和DIPIN结果，生成标准化报告"""
    input:
        sma_workout = rules.sma_analysis.output.sma_workout,
        dipin_alpha = rules.dipin_analysis.output.dipin_alpha,
        dipin_beta = rules.dipin_analysis.output.dipin_beta,
        dipin_other = rules.dipin_analysis.output.dipin_other,
        dirs_ready = rules.create_directories.output.directory_created
    output:
        sma_csv = "results/{sed_id}/{sample}/sma/{sample}_SMA_result.csv",
        dipin_csv = "results/{sed_id}/{sample}/dipin/{sample}_dipin_result.csv"
    params:
        workdir = os.path.abspath("."),
        uid = os.getuid(),
        gid = os.getgid(),
        race = lambda wildcards: SAMPLES.get(wildcards.sample, {}).get('Race', 'Asian'),  # 直接从样本配置获取种族信息，默认Asian
        docker_image = config['docker']['image'],
        docker_volumes = config['docker']['volumes'],
        script_path = config['Docker_scripts']['deal_sma_dipin_result']
    log:
        "results/{sed_id}/{sample}/logs/deal_sma_dipin_result.log"
    conda: "../env_yaml/xiaohaima_x03.yaml"
    shell:
        """
        echo "开始处理样本 {wildcards.sample} 的SMA和DIPIN结果..."
                
        # 运行结果处理脚本，直接传递种族信息
        docker run --rm \
            --user {params.uid}:{params.gid} \
            -v {params.workdir}:/workspace \
            {params.docker_volumes} \
            {params.docker_image} \
            python /scripts/deal_sma_dipin_result.py all \
                --alpha {input.dipin_alpha} \
                --beta {input.dipin_beta} \
                --other {input.dipin_other} \
                --sma {input.sma_workout} \
                --sample {wildcards.sample} \
                --output results/{wildcards.sed_id}/{wildcards.sample} \
                --race {params.race}
        
        
        echo "样本{wildcards.sample}————SMA和DIPIN结果处理完成"
        echo "样本种族信息: {params.race}"
        """
# ==================== 上传miniwes结果文件到海普云 ====================

rule upload_miniwes_results:
    """上传miniwes分析结果文件到海普云"""
    input:
        sma_csv = rules.deal_sma_dipin_result.output.sma_csv,
        dipin_csv = rules.deal_sma_dipin_result.output.dipin_csv
    output:
        upload_success = "results/{sed_id}/{sample}/upload_miniwes_success.txt"
    params:
        # 获取前6个字符，为批次的年月
        sed_id_short = lambda wildcards: wildcards.sed_id[:6],
        ossutil_path = config['upload']['ossutil_path'],
        oss_bucket =config['upload']['oss_bucket']
    log:
        "results/{sed_id}/{sample}/logs/upload_miniwes_results.log"
    conda: "../env_yaml/xiaohaima_x03.yaml"
    shell:
        """
        echo "开始上传样本 {wildcards.sample} 的miniwes分析结果文件..."
        
        # 上传SMA结果文件
        {params.ossutil_path} cp -f {input.sma_csv} {params.oss_bucket}/{params.sed_id_short}/{wildcards.sed_id}/
        
        # 上传DIPIN结果文件
        {params.ossutil_path} cp -f {input.dipin_csv} {params.oss_bucket}/{params.sed_id_short}/{wildcards.sed_id}/
        
        # 创建上传成功标志文件
        echo "MinIWES upload completed at $(date)" > {output.upload_success}
        
        echo "样本{wildcards.sample}————miniwes结果文件上传成功"
        """