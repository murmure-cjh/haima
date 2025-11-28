# core_analysis.smk - 核心分析模块（fastp、比对、去重、重比对、变异检测等）
# 功能：执行cfDNA测序数据的完整分析流程，包括质量控制、比对、去重、变异检测等
# 依赖：需要配置config.yaml文件和相关环境

# ==================== 配置和包含文件 ====================
configfile: "config.yaml"

# 包含通用函数定义
include: "common_def.smk"

# ==================== 目录创建规则 ====================
rule create_directories:
    """创建样本专用的目录结构"""
    output:
        directory_created = touch("results/{sed_id}/{sample}/.directory_created")
    params:
        # 使用lambda函数动态生成目录创建命令
        create_cmd = lambda wildcards: create_sample_dirs_shell(wildcards)
    log:
        "results/{sed_id}/{sample}/logs/create_directories.log"
    conda: "../env_yaml/xiaohaima_x03.yaml"
    shell:
        """
        echo "开始创建样本 {wildcards.sample} 的目录结构..."
        {params.create_cmd}
        # 创建标记文件表示目录结构已就绪
        touch {output.directory_created}
        echo "目录结构创建完成"
        """

# ==================== FASTP质控规则 ====================
rule fastp:
    """使用fastp进行原始测序数据的质量控制和过滤"""
    input:
        fq1 = lambda wildcards: get_fastq_path(wildcards, 'R1'),
        fq2 = lambda wildcards: get_fastq_path(wildcards, 'R2'),
        dirs_ready = rules.create_directories.output.directory_created
    output:
        clean_fq1 = "results/{sed_id}/{sample}/qc/{sample}_R1.fastq.gz",
        clean_fq2 = "results/{sed_id}/{sample}/qc/{sample}_R2.fastq.gz",
        html = "results/{sed_id}/{sample}/qc/{sample}_cfDNA_fastp.html",
        json = "results/{sed_id}/{sample}/qc/{sample}_cfDNA_fastp.json"
    params:
        # fastp质量控制参数：
        fastp_params = "-q 15 -u 40 -3 -W 4 -M 20 -n 5 -l 15 -w 10",
        workdir = os.path.abspath("."),
        uid = os.getuid(),
        gid = os.getgid(),
        docker_image = config['docker']['image'],
        docker_volumes = config['docker']['volumes']
    log:
        "results/{sed_id}/{sample}/logs/qc_fastp.log"
    threads: config['threads']
    conda: "../env_yaml/xiaohaima_x03.yaml"
    shell:
        """
        echo "开始对样本 {wildcards.sample} 进行fastp质控..."
        
        # 使用Docker容器运行fastp
        docker run --rm \
            --user {params.uid}:{params.gid} \
            -v {params.workdir}:/workspace \
            {params.docker_volumes} \
            {params.docker_image} \
            fastp {params.fastp_params} \
                -i {input.fq1} -I {input.fq2} \
                -o {output.clean_fq1} -O {output.clean_fq2} \
                -h {output.html} -j {output.json} \
                2>{log}
        
        echo "样本{wildcards.sample}————fastp质控完成，输出文件: {output.clean_fq1}, {output.clean_fq2}"
        """

# ==================== SENTIEON BWA比对规则 ====================
rule bwa_mem:
    """使用Sentieon优化的BWA MEM进行序列比对和排序"""
    input:
        fq1 = rules.fastp.output.clean_fq1,
        fq2 = rules.fastp.output.clean_fq2,
        dirs_ready = rules.create_directories.output.directory_created
    output:
        bam = "results/{sed_id}/{sample}/alignment/{sample}.sort.bam"
    params:
        sentieon_dir = config['sentieon_dir'],
        sentieon_license = config['sentieon_license'],
        ref = config['ref']
    log:
        "results/{sed_id}/{sample}/logs/bwa.log"
    threads: config['threads']
    conda: "../env_yaml/xiaohaima_x03.yaml"
    shell:
        """
        echo "开始对样本 {wildcards.sample} 进行BWA比对..."
        
        # 设置Sentieon许可证环境变量
        export SENTIEON_LICENSE={params.sentieon_license}
        
        # 执行BWA比对和排序流程
        {params.sentieon_dir} bwa mem \
            -R "@RG\\tID:{wildcards.sample}\\tLB:{wildcards.sample}\\tPL:ILLUMINA\\tSM:{wildcards.sample}" \
            -t {threads} -k 32 -M {params.ref} \
            {input.fq1} {input.fq2} | \
        {params.sentieon_dir} util sort \
            -o {output.bam} -t {threads} --sam2bam -i -
        
        echo "样本{wildcards.sample}————BWA比对完成，输出文件: {output.bam}"
        """

# ==================== 去重规则 ====================
rule dedup:
    """使用Sentieon进行PCR重复标记和去除"""
    input:
        bam = rules.bwa_mem.output.bam,
        dirs_ready = rules.create_directories.output.directory_created
    output:
        dedup_bam = "results/{sed_id}/{sample}/alignment/{sample}.dedup.bam",
        score = "results/{sed_id}/{sample}/tmp/{sample}.score.txt",
        metrics = "results/{sed_id}/{sample}/alignment/{sample}.dedup_metrics.txt"
    params:
        sentieon_dir = config['sentieon_dir'],
        sentieon_license = config['sentieon_license'],
        ref = config['ref']
    log:
        "results/{sed_id}/{sample}/logs/dedup.log"
    threads: config['threads']
    conda: "../env_yaml/xiaohaima_x03.yaml"
    shell:
        """
        echo "开始对样本 {wildcards.sample} 进行去重处理..."
        
        export SENTIEON_LICENSE={params.sentieon_license}
        
        # 步骤1: 计算位点收集器得分信息
        {params.sentieon_dir} driver \
            -r {params.ref} -t {threads} -i {input.bam} \
            --algo LocusCollector --fun score_info {output.score}
        
        # 步骤2: 执行去重算法
        {params.sentieon_dir} driver \
            -t {threads} -i {input.bam} \
            --algo Dedup --rmdup --score_info {output.score} \
            --metrics {output.metrics} {output.dedup_bam} \
            2>{output.dedup_bam}.dup.log
        
        echo "样本{wildcards.sample}————去重处理完成，输出文件: {output.dedup_bam}"
        """

# ==================== 重比对和BQSR规则 ====================
rule realign_bqsr:
    """执行局部重比对和碱基质量值重校准"""
    input:
        dedup_bam = rules.dedup.output.dedup_bam,
        dirs_ready = rules.create_directories.output.directory_created
    output:
        realigned_bam = "results/{sed_id}/{sample}/alignment/{sample}.realigned.bam",
        bqsr_bam = "results/{sed_id}/{sample}/alignment/{sample}.bqsr.bam",
        recal_table = "results/{sed_id}/{sample}/alignment/{sample}.recal_data.table"
    params:
        sentieon_dir = config['sentieon_dir'],
        sentieon_license = config['sentieon_license'],
        ref = config['ref'],
        mills_indel = config['mills_indel'],
        phase1_indel = config['phase1_indel'],
        bed = get_bed_file,  # 使用动态BED文件获取函数
        dipin_bed_slop60 = config['dipin_bed_slop60']
    log:
        realign = "results/{sed_id}/{sample}/logs/{sample}_realign.log",
        haplotyper = "results/{sed_id}/{sample}/logs/{sample}_Haplotyper.log"
    threads: config['threads']
    conda: "../env_yaml/xiaohaima_x03.yaml"
    shell:
        """
        echo "开始对样本 {wildcards.sample} 进行重比对和BQSR..."
        
        export SENTIEON_LICENSE={params.sentieon_license}
        
        # 步骤1: 局部重比对（针对indel区域）
        {params.sentieon_dir} driver \
            -r {params.ref} -t {threads} -i {input.dedup_bam} \
            --algo Realigner \
            -k {params.mills_indel} -k {params.phase1_indel} \
            --interval_list {params.bed} \
            {output.realigned_bam} \
            2> {log.realign}
        
        # 步骤2: 碱基质量值重校准（BQSR）
        {params.sentieon_dir} driver \
            -t {threads} -r {params.ref} -i {output.realigned_bam} \
            --interval {params.bed} \
            --algo QualCal \
            -k {params.mills_indel} -k {params.phase1_indel} \
            {output.recal_table}
        
        # 步骤3: 应用BQSR校正
        {params.sentieon_dir} driver \
            -t {threads} -r {params.ref} -i {output.realigned_bam} \
            --read_filter QualCalFilter,table={output.recal_table},indel=false \
            --interval {params.bed} \
            --algo ReadWriter \
            {output.bqsr_bam}
        
        echo "样本{wildcards.sample}————重比对和BQSR完成，输出文件: {output.bqsr_bam}"
        """

# ==================== 覆盖度分析规则 ====================
rule bamdst:
    """使用bamdst进行目标区域覆盖度分析"""
    input:
        dedup_bam = rules.dedup.output.dedup_bam,
        fastp_json = rules.fastp.output.json,
        dirs_ready = rules.create_directories.output.directory_created
    output:
        panel_coverage = "results/{sed_id}/{sample}/qc/coverage.report",
        panel_depth_plot = "results/{sed_id}/{sample}/qc/depth_distribution.plot",
    params:
        bamdst_cmd = config['bamdst'],
        cutoffdepth = "[1,5,10,20,30,50,100]", 
        panelfile = get_bed_qc_file, # 面板区域BED文件
        # targetfile = config['haima_bed'], # 目标区域BED文件
        workdir = os.path.abspath("."),
        uid = os.getuid(),
        gid = os.getgid(),
        docker_image = config['docker']['image'],
        docker_volumes = config['docker']['volumes']
    log:
        "results/{sed_id}/{sample}/logs/bamdst.log"
    threads: config['threads']
    conda: "../env_yaml/xiaohaima_x03.yaml"
    shell:
        """
        echo "开始对样本 {wildcards.sample} 进行覆盖度分析..."
        
        # 分析面板区域的覆盖度
        docker run --rm \
            --user {params.uid}:{params.gid} \
            -v {params.workdir}:/workspace \
            {params.docker_volumes} \
            {params.docker_image} \
            bamdst --cutoffdepth {params.cutoffdepth} \
                -p {params.panelfile} \
                -o results/{wildcards.sed_id}/{wildcards.sample}/qc \
                {input.dedup_bam}
               
        echo "样本{wildcards.sample}————覆盖度分析完成"
        """

# ==================== 质控分析规则 ====================
# 因为最终质控报告的文件结构问题，会使用到target_coverage，后续如果修改，需要和IT解读沟通之后，再修改流程/脚本和删除对应的配置 
rule qc_analysis:
    """整合各项质控指标生成样本质控报告"""
    input:
        dedup_bam = rules.dedup.output.dedup_bam,
        dedup_metrics = rules.dedup.output.metrics,
        fastp_json = rules.fastp.output.json,
        panel_coverage = rules.bamdst.output.panel_coverage,
        panel_depth_plot = rules.bamdst.output.panel_depth_plot,
        dirs_ready = rules.create_directories.output.directory_created
    output:
        sample_info_json = "results/{sed_id}/{sample}/qc/{sample}.sample_information.json",
        qc_info_txt = "results/{sed_id}/{sample}/qc/{sample}.sample_information.txt",
        qc_info_csv = "results/{sed_id}/{sample}/qc/{sample}_QC_result.csv",
    params:
        panelfile = get_bed_qc_file,
        targetfile = config['haima_bed'],
        sampleid = get_sample_id,
        workdir = os.path.abspath("."),
        uid = os.getuid(),
        gid = os.getgid(),
        docker_image = config['docker']['image'],
        docker_volumes = config['docker']['volumes'],
        script_path = config['Docker_scripts']['qc_fload_analysis']
    log:
        "results/{sed_id}/{sample}/logs/qc_analysis.log"
    conda: "../env_yaml/xiaohaima_x03.yaml"
    shell:
        """
        echo "开始生成样本 {wildcards.sample} 的质控报告..."
        
        # 运行Python质控分析脚本
        docker run --rm \
            --user {params.uid}:{params.gid} \
            -v {params.workdir}:/workspace \
            {params.docker_volumes} \
            {params.docker_image} \
            python {params.script_path} \
                {wildcards.sample} \
                {input.dedup_bam} \
                {input.dedup_metrics} \
                {params.panelfile} \
                results/{wildcards.sed_id}/{wildcards.sample}/qc \
                {params.sampleid} \
                2>{log}
        
        echo "样本{wildcards.sample}————质控报告生成完成"
        """

# ==================== 变异检测规则 ====================
rule haplotyper:
    """使用Sentieon Haplotyper进行变异检测"""
    input:
        bqsr_bam = rules.realign_bqsr.output.bqsr_bam,
        dirs_ready = rules.create_directories.output.directory_created
    output:
        raw_vcf = "results/{sed_id}/{sample}/variant/{sample}.raw.vcf"
    params:
        sentieon_dir = config['sentieon_dir'],
        sentieon_license = config['sentieon_license'],
        ref = config['ref'],
        bed = get_bed_file,  # 使用动态BED文件
        dipin_bed_slop60 = config['dipin_bed_slop60']
    log:
        "results/{sed_id}/{sample}/logs/{sample}_Haplotyper.log"
    threads: config['threads']
    conda: "../env_yaml/xiaohaima_x03.yaml"
    shell:
        """
        echo "开始对样本 {wildcards.sample} 进行变异检测..."
        
        export SENTIEON_LICENSE={params.sentieon_license}
        
        # 执行Haplotyper变异检测
        {params.sentieon_dir} driver \
            -r {params.ref} -t {threads} -i {input.bqsr_bam} \
            --interval {params.bed} \
            --algo Haplotyper \
            --emit_conf 30.0 --call_conf 10 \
            {output.raw_vcf} \
            2>> {log}
        
        echo "样本{wildcards.sample}————变异检测完成，输出文件: {output.raw_vcf}"
        """

# ==================== 性别检测规则 ====================
# 小天机因为不包含性染色体，所以不进行分析
rule sample_gender_analysis:
    """基于染色体覆盖度进行样本性别检测"""
    input:
        dedup_bam = rules.dedup.output.dedup_bam,
        dirs_ready = rules.create_directories.output.directory_created
    output:
        gender_file = "results/{sed_id}/{sample}/gender/{sample}.gender.txt"
    params:
        panel = get_panel_type,
        gender_params = lambda wildcards: get_gender_params(wildcards),
        workdir = os.path.abspath("."),
        uid = os.getuid(),
        gid = os.getgid(),
        docker_image = config['docker']['image'],
        docker_volumes = config['docker']['volumes']
    log:
        "results/{sed_id}/{sample}/logs/sample_gender_analysis.log"
    conda: "../env_yaml/xiaohaima_x03.yaml"
    shell:
        """
        echo "开始进行样本 {wildcards.sample} 的性别检测..."
        
        panel="{params.panel}"
        
        # 检查panel类型，NAD-XTJ2211类型跳过性别分析
        if [ "$panel" = "NAD-XTJ2211" ]; then
            echo "跳过性别分析，panel类型为NAD-XTJ2211"
            echo "Unknown" > {output.gender_file}
        else
            # 运行ngs-bits SampleGender分析
            docker run --rm \
                --user {params.uid}:{params.gid} \
                -v {params.workdir}:/workspace \
                {params.docker_volumes} \
                {params.docker_image} \
                SampleGender -in {input.dedup_bam} -method xy {params.gender_params} > {output.gender_file}
        fi
        
        echo "样本{wildcards.sample}————性别检测完成，结果文件: {output.gender_file}"
        """
 
# ==================== 变异过滤规则 ====================
rule filter_variants:
    """对原始变异进行质量过滤和分类"""
    input:
        raw_vcf = rules.haplotyper.output.raw_vcf,
        dirs_ready = rules.create_directories.output.directory_created
    output:
        snp_vcf = "results/{sed_id}/{sample}/variant/{sample}.raw.snp.vcf",
        indel_vcf = "results/{sed_id}/{sample}/variant/{sample}.raw.indel.vcf",
        multi_vcf = "results/{sed_id}/{sample}/variant/{sample}.raw.multi.vcf",
        filtered_vcf = "results/{sed_id}/{sample}/variant/{sample}.filter.vcf"
    params:
        ref = config['ref'],
        workdir = os.path.abspath("."),
        uid = os.getuid(),
        gid = os.getgid(),
        docker_image = config['docker']['image'],
        docker_volumes = config['docker']['volumes'],
        # SNP过滤条件
        snp_filter = '"QD >= 2.0; FS <= 60.0; MQ >= 40.0; SOR <= 4.0; MQRankSum >= -12.5; ReadPosRankSum >= -8.0"',
        # Indel过滤条件
        indel_filter = '"QD >= 2.0; FS <= 200.0; ReadPosRankSum >= -20.0; SOR <= 10.0"',
        # 多等位基因过滤条件，这部分的过滤条件有优化的空间
        multi_filter = '"QD >= 5.0; FS <= 30.0; MQ >= 50.0; SOR <= 3.0; MQRankSum >= -8.0; ReadPosRankSum >= -5.0"' 
    log:
        "results/{sed_id}/{sample}/logs/filter_variants.log"
    conda: "../env_yaml/xiaohaima_x03.yaml"
    shell:
        """
        echo "开始对样本 {wildcards.sample} 的变异进行过滤..."
        
        # 步骤1: 分离SNP变异
        docker run --rm \
            --user {params.uid}:{params.gid} \
            -v {params.workdir}:/workspace \
            {params.docker_volumes} \
            {params.docker_image} \
            VcfFilter -ref {params.ref} -in {input.raw_vcf} -variant_type snp -out {output.snp_vcf}
        
        # 步骤2: 分离Indel变异
        docker run --rm \
            --user {params.uid}:{params.gid} \
            -v {params.workdir}:/workspace \
            {params.docker_volumes} \
            {params.docker_image} \
            VcfFilter -ref {params.ref} -in {input.raw_vcf} -variant_type indel -out {output.indel_vcf}
        
        # 步骤3: 分离多等位基因变异
        docker run --rm \
            --user {params.uid}:{params.gid} \
            -v {params.workdir}:/workspace \
            {params.docker_volumes} \
            {params.docker_image} \
            VcfFilter -ref {params.ref} -in {input.raw_vcf} -variant_type multi-allelic -out {output.multi_vcf}
        
        # 步骤4: 对SNP进行质量过滤
        docker run --rm \
            --user {params.uid}:{params.gid} \
            -v {params.workdir}:/workspace \
            {params.docker_volumes} \
            {params.docker_image} \
            VcfFilter -ref {params.ref} -in {output.snp_vcf} -out {output.snp_vcf}.tmp \
                -info {params.snp_filter} \
                -filter_empty
        
        # 步骤5: 对Indel进行质量过滤
        docker run --rm \
            --user {params.uid}:{params.gid} \
            -v {params.workdir}:/workspace \
            {params.docker_volumes} \
            {params.docker_image} \
            VcfFilter -ref {params.ref} -in {output.indel_vcf} -out {output.indel_vcf}.tmp \
                -info {params.indel_filter} \
                -filter_empty
        
        # 步骤6: 对多等位基因进行质量过滤
        docker run --rm \
            --user {params.uid}:{params.gid} \
            -v {params.workdir}:/workspace \
            {params.docker_volumes} \
            {params.docker_image} \
            VcfFilter -ref {params.ref} -in {output.multi_vcf} -out {output.multi_vcf}.tmp \
                -info {params.multi_filter} \
                -filter_empty
        
        # 步骤7: 合并所有过滤后的变异
        # 7.1 写入VCF头信息（从SNP文件中获取）
        grep '^#' {output.snp_vcf}.tmp > {output.filtered_vcf}
        # 7.2 合并所有通过过滤的变异行
        cat {output.snp_vcf}.tmp {output.indel_vcf}.tmp {output.multi_vcf}.tmp | grep -v '^#' >> {output.filtered_vcf}
        
        # 步骤8: 清理临时文件
        rm {output.snp_vcf}.tmp {output.indel_vcf}.tmp {output.multi_vcf}.tmp
        
        echo "样本{wildcards.sample}————变异过滤完成，输出文件: {output.filtered_vcf}"
        """