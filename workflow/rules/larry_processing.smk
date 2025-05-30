rule extract_barcodes:
    input: 
        fb = "data/lane_merged/{{sample}}_FB_S1_L001_{}_001.fastq.gz".format(LARRY["read_feature_bc"]),
        cb = "data/lane_merged/{{sample}}_FB_S1_L001_{}_001.fastq.gz".format(LARRY["read_cellular_bc"]),
    output:
        filt_fb = temp(
            expand(
                "data/bc_filt/{{sample}}_FB_S1_L001_{read_fb}_001_{larry_color}.fastq.gz", 
                larry_color = LARRY_COLORS, 
                read_fb = LARRY["read_feature_bc"]
                )
            ),
        filt_cb = temp("data/clean/{{sample}}_FB_S1_L001_{}_001.fastq.gz".format(LARRY["read_cellular_bc"]))
    params:
        barcode_dict = LARRY["bc_patterns"]
    log:
        "results/00_logs/extract_barcodes/{sample}.log"
    benchmark:
        "results/benchmarks/extract_barcodes/{sample}.txt"
    container:
         config["singularity"]["seurat_sif"]
    resources:
        mem_mb  = get_mem_mb(RESOURCES["extract_barcodes"]["mem_mb"], 20000),
        runtime = RESOURCES["extract_barcodes"]["runtime"]
    retries:
        RESOURCES["extract_barcodes"]["retries"]
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
        mem_mb  = get_mem_mb(RESOURCES["collapse_fastq_hd"]["mem_mb"], 20000),
        runtime = RESOURCES["collapse_fastq_hd"]["runtime"]
    retries:
        RESOURCES["collapse_fastq_hd"]["retries"]
    container:
        config["singularity"]["umicollapse"]
    shell:
        """
        java -Xms{resources.mem_mb}m -Xmx{resources.mem_mb}m -Xss1G -jar /UMICollapse/umicollapse.jar fastq -k {wildcards.hd} --tag -i {input} -o {output} 2> {log}
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
        mem_mb  = get_mem_mb(RESOURCES["correct_barcodes"]["mem_mb"], 10000),
        runtime = RESOURCES["correct_barcodes"]["runtime"]
    retries:
        RESOURCES["correct_barcodes"]["retries"]
    container:
         config["singularity"]["seurat_sif"]
    script:
        "../scripts/python/correct_barcodes.py"


rule merge_corrected_fastq:
    input:
        expand(
            "data/corrected/{{sample}}_FB_S1_L001_{{read_fb}}_001_{larry_color}_collapsed-hd{hd}-corrected.fastq.gz",
            larry_color = LARRY_COLORS,
            hd          = LARRY["hamming_distance"]
        )
    output:
        temp("data/clean/{sample}_FB_S1_L001_{read_fb}_001.fastq.gz")
    log:
        "results/00_logs/merge_corrected_fastq/{sample}_{read_fb}.log"
    benchmark:
        "results/benchmarks/merge_corrected_fastq/{sample}_{read_fb}.txt"
    shell:
        """
        cat {input} > {output}
        """


rule generate_feature_ref_larry:
    input:
        expand(
            "data/feature_reference/{sample}_{read_fb}_{larry_color}_hd{hd}_feature_reference.csv", 
            sample      = SAMPLES,
            read_fb     = LARRY["read_feature_bc"],
            larry_color = LARRY_COLORS,
            hd          = LARRY["hamming_distance"]
        )
    output:
        "data/feature_reference/Feature_reference_larry.csv"
    log:
        "results/00_logs/generate_feature_ref/log"
    container:
        "../envs/Seurat.yaml"
    script:
        "../scripts/R/generate_feature_ref_larry.R"