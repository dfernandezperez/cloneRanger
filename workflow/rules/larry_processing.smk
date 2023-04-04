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

rule extract_barcodes:
    input: 
        "data/{sample}_FB_S1_L001_{read}_001.fastq.gz"
    output:
        parsed_fq   = "data/bc_filt/{sample}_FB_S1_L001_{read}_001.fastq.gz",
        feature_ref = "results/00_logs/feature_bc_libraries/{sample}_{read}_feature_reference.csv"
    params:
        barcode_dict    = config["feature_bc_config"]["bc_patterns"],
        read_feature_bc = config["feature_bc_config"]["read_feature_bc"]
    log:
        "results/00_logs/extract_barcodes/{sample}_{read}.log"
    benchmark:
        "results/benchmarks/extract_barcodes/{sample}_{read}.txt"
    script:
        "../scripts/extract_barcodes.py"

rule extract_barcodes_v2:
    input: 
        "data/{sample}_FB_S1_L001_{read}_001.fastq.gz"
    output:
        parsed_fq = expand("data/bc_filt/{{sample}}_FB_S1_L001_{{read}}_001_{larry_color}.fastq.gz", larry_color = LARRY_COLORS)
    params:
        barcode_dict    = config["feature_bc_config"]["bc_patterns"],
        read_feature_bc = config["feature_bc_config"]["read_feature_bc"]
    log:
        "results/00_logs/extract_barcodes/{sample}_{read}.log"
    benchmark:
        "results/benchmarks/extract_barcodes/{sample}_{read}.txt"
    script:
        "../scripts/extract_barcodes_v2.py"

rule collapse_fastq_hd:
    input:
        expand("data/bc_filt/{sample}_FB_S1_L001_{read}_001.fastq.gz", read = config["feature_bc_config"]["read_feature_bc"], sample = SAMPLES)
    output:
        "data/collapsed_barcodes_hd{hd}.fastq.gz"
    log:
        "results/00_logs/collapse_fastq_hd{hd}/log"
    benchmark:
        "results/benchmarks/hd_collapse/hamming_distance_{hd}.txt"
    shell:
        """
        cat {input} > tmp.fq.gz
        /home/dfernandezp/miniconda3/bin/java -Xmx200G -Xss1G -jar /stemcell/scratch/dfernandezp/IndranilSingh/drug_screening_rabseq/software/UMICollapse/umicollapse.jar fastq -k {wildcards.hd} --tag -i tmp.fq.gz -o {output}
        rm tmp.fq.gz
        """

rule generate_feature_ref:
    input:
        expand("results/00_logs/feature_bc_libraries/{sample}_feature_reference.csv", sample = SAMPLES)
    output:
        "results/00_logs/feature_bc_libraries/Feature_reference.csv"
    log:
        "results/00_logs/generate_feature_ref/log"
    script:
        "../scripts/generate_feature_ref.R"