rule extract_barcodes:
    input: 
        fb = "data/lane_merged/{{sample}}_FB_S1_L001_{}_001.fastq.gz".format(config["feature_bc_config"]["read_feature_bc"]),
        cb = "data/lane_merged/{{sample}}_FB_S1_L001_{}_001.fastq.gz".format(config["feature_bc_config"]["read_cellular_bc"]),
    output:
        filt_fb = temp(
            expand(
                "data/bc_filt/{{sample}}_FB_S1_L001_{read_fb}_001_{larry_color}.fastq.gz", 
                larry_color = LARRY_COLORS, 
                read_fb = config["feature_bc_config"]["read_feature_bc"]
                )
            ),
        filt_cb = "data/clean/{{sample}}_FB_S1_L001_{}_001.fastq.gz".format(config["feature_bc_config"]["read_cellular_bc"])
    params:
        barcode_dict = config["feature_bc_config"]["bc_patterns"]
    log:
        "results/00_logs/extract_barcodes/{sample}.log"
    benchmark:
        "results/benchmarks/extract_barcodes/{sample}.txt"
    conda:
         "../envs/python.yaml"
    resources:
        mem_mb = RESOURCES["extract_barcodes"]["MaxMem"]
    script:
        "../scripts/python/extract_barcodes.py"


rule collapse_fastq_hd:
    input:
        "data/bc_filt/{sample}_FB_S1_L001_{read_fb}_001_{larry_color}.fastq.gz"
    output:
        temp("data/collapsed/{sample}_FB_S1_L001_{read_fb}_001_{larry_color}_collapsed-hd{hd}.fastq.gz")
    log:
        "results/00_logs/collapse_fastq_hd/{sample}_{read_fb}_{larry_color}_{hd}.log"
    benchmark:
        "results/benchmarks/collapse_fastq_hd/{sample}_{read_fb}_{larry_color}_{hd}.txt"
    resources:
        mem_mb = RESOURCES["collapse_fastq_hd"]["MaxMem"]
    container:
        "docker://dfernand/umicollapse:07506e661496033ca059e04a7771633b7b70721f"
    shell:
        """
        java -Xmx200G -Xss1G -jar /UMICollapse/umicollapse.jar fastq -k {wildcards.hd} --tag -i {input} -o {output} 2> {log}
        """


rule correct_barcodes:
    input:
        "data/collapsed/{sample}_FB_S1_L001_{read_fb}_001_{larry_color}_collapsed-hd{hd}.fastq.gz"
    output:
        corrected_fq = temp("data/corrected/{sample}_FB_S1_L001_{read_fb}_001_{larry_color}_collapsed-hd{hd}-corrected.fastq.gz"),
        feature_ref  = temp("data/feature_reference/{sample}_{read_fb}_{larry_color}_hd{hd}_feature_reference.csv")
    log:
        "results/00_logs/correct_barcodes/{sample}_{read_fb}_{larry_color}_{hd}.log"
    benchmark:
        "results/benchmarks/correct_barcodes/{sample}_{read_fb}_{larry_color}_{hd}.txt"
    resources:
        mem_mb = RESOURCES["correct_barcodes"]["MaxMem"]
    conda:
         "../envs/python.yaml"
    script:
        "../scripts/python/correct_barcodes.py"


rule merge_corrected_fastq:
    input:
        expand(
            "data/corrected/{{sample}}_FB_S1_L001_{{read_fb}}_001_{larry_color}_collapsed-hd{hd}-corrected.fastq.gz",
            larry_color = LARRY_COLORS,
            hd          = config["feature_bc_config"]["hamming_distance"]
        )
    output:
        "data/clean/{sample}_FB_S1_L001_{read_fb}_001.fastq.gz"
    log:
        "results/00_logs/merge_corrected_fastq/{sample}_{read_fb}.log"
    benchmark:
        "results/benchmarks/merge_corrected_fastq/{sample}_{read_fb}.txt"
    shell:
        """
        cat {input} > {output}
        """


rule generate_feature_ref:
    input:
        expand(
            "data/feature_reference/{sample}_{read_fb}_{larry_color}_hd{hd}_feature_reference.csv", 
            sample      = SAMPLES,
            read_fb     = config["feature_bc_config"]["read_feature_bc"],
            larry_color = LARRY_COLORS,
            hd          = config["feature_bc_config"]["hamming_distance"]
        )
    output:
        "data/feature_reference/Feature_reference.csv"
    log:
        "results/00_logs/generate_feature_ref/log"
    conda:
        "../envs/Seurat.yaml"
    script:
        "../scripts/R/generate_feature_ref.R"