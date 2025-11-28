# annotation_analysis.smk - 注释和分析流程模块
# 功能：对变异检测结果进行功能注释、SNP分析和结果格式转换
# 依赖：需要config.yaml配置文件和common_def.smk通用函数

# ==================== 配置和包含文件 ====================
configfile: "config.yaml"

# 包含通用函数定义
include: "common_def.smk"

# ==================== ANNOVAR和VEP注释规则 ====================
rule annovar_vep_annotation:
    """使用ANNOVAR和VEP对过滤后的变异进行功能注释"""
    input:
        vcf = rules.filter_variants.output.filtered_vcf,
        dirs_ready = rules.create_directories.output.directory_created,
        config_file = "config.yaml" 
    output:
        anno_txt = "results/{sed_id}/{sample}/annotation/{sample}_anno.txt",
        anno_json = "results/{sed_id}/{sample}/annotation/{sample}_anno.json"
    params:
        sample = "{sample}",
        outdir = get_sample_annotation_dir,
        workdir = os.path.abspath("."),
        uid = os.getuid(),
        gid = os.getgid(),
        docker_image = config['docker']['image'],
        docker_volumes = config['docker']['volumes'],
        script_path = config['Docker_scripts']['anno_caller']
    log:
        "results/{sed_id}/{sample}/logs/annovar_vep.log"
    conda: "../env_yaml/xiaohaima_x03.yaml"
    shell:
        """
        echo "开始对样本 {wildcards.sample} 进行ANNOVAR和VEP注释..."
        
        # 使用Docker容器运行注释脚本
        docker run --rm \
            --user {params.uid}:{params.gid} \
            -v {params.workdir}:/workspace \
            {params.docker_volumes} \
            {params.docker_image} \
            python {params.script_path} \
                --vcf {input.vcf} \
                --sample {params.sample} \
                --output {params.outdir} \
                --config {input.config_file} \
                2>{log}
        
        echo "样本{wildcards.sample}————ANNOVAR和VEP注释完成，输出文件: {output.anno_txt}, {output.anno_json}"
        """

# ==================== 海马结果处理规则 ====================
rule haimaresult_processing:
    """对注释结果进行后处理"""
    input:
        anno_json = rules.annovar_vep_annotation.output.anno_json,
        dirs_ready = rules.create_directories.output.directory_created,
        config_file = "config.yaml"  
    output:
        haima_json = "results/{sed_id}/{sample}/annotation/{sample}_haima.json"
    params:
        sample = "{sample}",
        workdir = os.path.abspath("."),
        uid = os.getuid(),
        gid = os.getgid(),
        docker_image = config['docker']['image'],
        docker_volumes = config['docker']['volumes'],
        script_path = config['Docker_scripts']['haimaresult']
    log:
        "results/{sed_id}/{sample}/logs/haimaresult.log"
    conda: "../env_yaml/xiaohaima_x03.yaml"
    shell:
        """
        echo "开始对样本 {wildcards.sample} 的注释结果进行海马格式处理..."
        
        # 运行海马结果处理脚本
        docker run --rm \
            --user {params.uid}:{params.gid} \
            -v {params.workdir}:/workspace \
            {params.docker_volumes} \
            {params.docker_image} \
            python {params.script_path} \
                -i {input.anno_json} \
                -o {output.haima_json} \
                --config {input.config_file} \
                2>{log}
        
        echo "样本{wildcards.sample}————海马结果处理完成，输出文件: {output.haima_json}"
        """

# ==================== SNP分析规则 ====================
rule snp_analysis:
    """使用haima_snp.py进行药物SNP相关分析"""
    input:
        bam = rules.realign_bqsr.output.bqsr_bam,
        anno_json = rules.annovar_vep_annotation.output.anno_json,
        dirs_ready = rules.create_directories.output.directory_created,
    output:
        snp_results = "results/{sed_id}/{sample}/snp/{sample}_snp_results.json"
    params:
        durg_bed = config['durg_bed'],  # 药物相关BED文件
        ref = config['ref'], 
        outdir = get_sample_snp_dir,   
        workdir = os.path.abspath("."),
        uid = os.getuid(),
        gid = os.getgid(),
        docker_image = config['docker']['image'],
        docker_volumes = config['docker']['volumes'],
        script_path = config['Docker_scripts']['haima_snp']
    log:
        "results/{sed_id}/{sample}/logs/snp_analysis.log"
    conda: "../env_yaml/xiaohaima_x03.yaml"
    shell:
        """
        echo "开始对样本 {wildcards.sample} 进行SNP分析..."
        
        # 运行SNP分析脚本
        docker run --rm \
            --user {params.uid}:{params.gid} \
            -v {params.workdir}:/workspace \
            {params.docker_volumes} \
            {params.docker_image} \
            python {params.script_path} \
                {wildcards.sample} \
                {input.bam} \
                {params.durg_bed} \
                {params.ref} \
                {params.outdir} \
                {input.anno_json} \
                2>{log}
        
        echo "样本{wildcards.sample}————SNP分析完成，输出文件: {output.snp_results}"
        """

# ==================== 海马结果转换规则 ====================
# 为匹配现有报告上传格式专门写的脚本
rule transform_haima_result:
    """转换海马结果格式以匹配原有系统格式要求"""
    input:
        haima_json = rules.haimaresult_processing.output.haima_json,
        snp_results = rules.snp_analysis.output.snp_results,
        dirs_ready = rules.create_directories.output.directory_created
    output:
        haima_germline_vcf = "results/{sed_id}/{sample}/annotation/{sample}_germline.vcf"
    params:
        sample = "{sample}",
        workdir = os.path.abspath("."),
        uid = os.getuid(),
        gid = os.getgid(),
        docker_image = config['docker']['image'],
        docker_volumes = config['docker']['volumes'],
        script_path = config['Docker_scripts']['transform_haima_result']
    log:
        "results/{sed_id}/{sample}/logs/transform_haima_result.log"
    conda: "../env_yaml/xiaohaima_x03.yaml"
    shell:
        """
        echo "开始转换样本 {wildcards.sample} 的海马结果格式..."
        
        # 运行格式转换脚本
        docker run --rm \
            --user {params.uid}:{params.gid} \
            -v {params.workdir}:/workspace \
            {params.docker_volumes} \
            {params.docker_image} \
            python {params.script_path} \
                -i {input.haima_json} \
                -a {input.snp_results} \
                -o {output.haima_germline_vcf} \
                2>{log}
        
        # 验证输出文件是否成功生成且非空
        echo "验证输出VCF文件..."
        if [ ! -s {output.haima_germline_vcf} ]; then
            echo "错误: 输出VCF文件为空或未生成" >&2
            exit 1
        fi
        
        echo "样本{wildcards.sample}————海马结果格式转换成功完成"
        echo "生成germline VCF文件: {output.haima_germline_vcf}"
        """

# ================= 上传结果文件到海普云 ====================
# 将质控结果和haima_germline_vcf结果上传到海普云
rule upload_file:
    """上传流程结果文件到海普云（非小全外步骤）"""
    input:
        haima_germline_vcf = rules.transform_haima_result.output.haima_germline_vcf,
        qc_info_txt = rules.qc_analysis.output.qc_info_txt
    output:
        # 添加上传成功的标志文件
        upload_success = "results/{sed_id}/{sample}/upload_success.txt"
    params:
        # 在这里获取前6个字符，为批次的年月
        sed_id_short = lambda wildcards: wildcards.sed_id[:6],
        ossutil_path = config['upload']['ossutil_path'],
        oss_bucket =config['upload']['oss_bucket']
    log:
        "results/{sed_id}/{sample}/logs/transform_haima_result.log"
    conda: "../env_yaml/xiaohaima_x03.yaml"
    shell:
        """
        echo "开始转换样本 {wildcards.sample} 开始上传文件..."
        # 上传文件
        {params.ossutil_path} cp -f {input.haima_germline_vcf} {params.oss_bucket}/{params.sed_id_short}/{wildcards.sed_id}/
        
        # 目前不上传质控信息
        #{params.ossutil_path} cp -f {input.qc_info_txt} {params.oss_bucket}/{params.sed_id_short}/{wildcards.sed_id}/
        
        # 创建上传成功标志文件
        echo "Upload completed at $(date)" > {output.upload_success}

        echo "样本{wildcards.sample}————结果文件上传成功"
        """