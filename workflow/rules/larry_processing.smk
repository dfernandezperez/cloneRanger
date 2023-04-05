rule create_library:
    output:
        "results/00_logs/feature_bc_libraries/{sample}_library.csv"
    params:
        feature_barcoding = get_library_type
    shell:
        """
        echo "fastqs,sample,library_type" > {output}
        echo "{params.feature_barcoding}" >> {output}
        """

# rule extract_barcodes:
#     input: 
#         "data/{sample}_FB_S1_L001_{read}_001.fastq.gz"
#     output:
#         parsed_fq   = "data/bc_filt/{sample}_FB_S1_L001_{read}_001.fastq.gz",
#         feature_ref = "results/00_logs/feature_bc_libraries/{sample}_{read}_feature_reference.csv"
#     params:
#         barcode_dict    = config["feature_bc_config"]["bc_patterns"],
#         read_feature_bc = config["feature_bc_config"]["read_feature_bc"]
#     log:
#         "results/00_logs/extract_barcodes/{sample}_{read}.log"
#     benchmark:
#         "results/benchmarks/extract_barcodes/{sample}_{read}.txt"
#     script:
#         "../scripts/extract_barcodes_old.py"

rule extract_barcodes:
    input: 
        "data/{sample}_FB_S1_L001_{read}_001.fastq.gz"
    output:
        parsed_fq = temp(expand("data/bc_filt/{{sample}}_FB_S1_L001_{{read}}_001_{larry_color}.fastq.gz", larry_color = LARRY_COLORS))
    params:
        barcode_dict    = config["feature_bc_config"]["bc_patterns"],
        read_feature_bc = config["feature_bc_config"]["read_feature_bc"]
    log:
        "results/00_logs/extract_barcodes/{sample}_{read}.log"
    benchmark:
        "results/benchmarks/extract_barcodes/{sample}_{read}.txt"
    script:
        "../scripts/extract_barcodes.py"

rule collapse_fastq_hd:
    input:
        "data/bc_filt/{sample}_FB_S1_L001_{read}_001_{larry_color}.fastq.gz"
    output:
        temp("data/bc_filt/collapsed/{sample}_FB_S1_L001_{read}_001_{larry_color}_collapsed-hd{hd}.fastq.gz")
    log:
        "results/00_logs/collapse_fastq_hd/{sample}_{read}_{larry_color}_{hd}.log"
    benchmark:
        "results/benchmarks/collapse_fastq_hd/{sample}_{read}_{larry_color}_{hd}.txt"
    shell:
        """
        /home/dfernandezp/miniconda3/bin/java -Xmx200G -Xss1G -jar /stemcell/scratch/dfernandezp/IndranilSingh/drug_screening_rabseq/software/UMICollapse/umicollapse.jar fastq -k {wildcards.hd} --tag -i {input} -o {output} 2> {log}
        """

rule correct_barcodes:
    input:
        "data/bc_filt/collapsed/{sample}_FB_S1_L001_{read}_001_{larry_color}_collapsed-hd{hd}.fastq.gz"
    output:
        corrected_fq = "data/bc_filt/collapsed/{sample}_FB_S1_L001_{read}_001_{larry_color}_collapsed-hd{hd}-corrected.fastq.gz",
        feature_ref  = "data/feature_reference/{sample}_{read}_{larry_color}_hd{hd}_feature_reference.csv"
    log:
        "results/00_logs/correct_barcodes/{sample}_{read}_{larry_color}_{hd}.log"
    benchmark:
        "results/benchmarks/correct_barcodes/{sample}_{read}_{larry_color}_{hd}.txt"
    script:
        "../scripts/correct_barcodes.py"

rule generate_feature_ref:
    input:
        expand(
            "data/feature_reference/{sample}_{read}_{larry_color}_hd{hd}_feature_reference.csv", 
            sample      = SAMPLES,
            read        = config["feature_bc_config"]["read_feature_bc"],
            larry_color = LARRY_COLORS,
            hd          = config["feature_bc_config"]["hamming_distance"]
        )
    output:
        "data/feature_reference/Feature_reference.csv"
    log:
        "results/00_logs/generate_feature_ref/log"
    script:
        "../scripts/generate_feature_ref.R"