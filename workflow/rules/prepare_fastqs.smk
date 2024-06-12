# This currently expects names as I receive them from the core
rule clean_names:
    input:
        get_fastqs,
    output:
        fw = temp("data/symlink/{sample}_{lib_type}_{lane}_R1.fastq.gz"),
        rv = temp("data/symlink/{sample}_{lib_type}_{lane}_R2.fastq.gz"),
    script:
        "../scripts/python/create_symlink.py"

rule clean_names_atac:
    input:
        get_atac_fastqs,
    output:
        temp("data/symlink/{sample}_{lib_type}_{lane}_R3.fastq.gz"),
    script:
        "../scripts/python/create_symlink_atac.py"

rule merge_lanes:
    input:
        fw = lambda w: expand(
            "data/symlink/{sample.sample_id}_{sample.lib_type}_{sample.lane}_R1.fastq.gz", 
            sample=units.loc[(w.sample, w.lib_type)].itertuples()
            ),
        rv = lambda w: expand(
            "data/symlink/{sample.sample_id}_{sample.lib_type}_{sample.lane}_R2.fastq.gz", 
            sample=units.loc[(w.sample, w.lib_type)].itertuples()
            )
    output:
        fw = temp("data/lane_merged/{sample}_{lib_type}_S1_L001_R1_001.fastq.gz"),
        rv = temp("data/lane_merged/{sample}_{lib_type}_S1_L001_R2_001.fastq.gz")
    log:
        "results/00_logs/merge_lanes/{sample}_{lib_type}.log"
    resources:
        runtime = RESOURCES["merge_lanes"]["runtime"]
    container:
        None
    shell:
        """
        cat {input.fw} > {output.fw} 2> {log}
        cat {input.rv} > {output.rv} 2>> {log}
        """

rule merge_lanes_atac:
    input:
        lambda w: expand(
            "data/symlink/{sample.sample_id}_{sample.lib_type}_S1_L00{sample.lane}_R3_001.fastq.gz", 
            sample=units.loc[(w.sample, w.lib_type)].itertuples()
            )
    output:
        temp("data/lane_merged/{sample}_{lib_type}_S1_L001_R3_001.fastq.gz")
    log:
        "results/00_logs/merge_lanes_atac/{sample}_{lib_type}.log"
    resources:
        runtime = RESOURCES["merge_lanes"]["runtime"]
    container:
        None
    shell:
        """
        cat {input} > {output} 2> {log}
        """
    
rule move_gex_fq:
    # Dummy rule to move the gex fastq files to the "clean" folder, in which the collapsed feature barcoding
    # fastq files will be stored.
    input:
        fw = "data/lane_merged/{sample}_GEX_S1_L001_R1_001.fastq.gz",
        rv = "data/lane_merged/{sample}_GEX_S1_L001_R2_001.fastq.gz"
    output:
        fw = temp("data/clean/{sample}_GEX_S1_L001_R1_001.fastq.gz"),
        rv = temp("data/clean/{sample}_GEX_S1_L001_R2_001.fastq.gz")
    container:
        None
    shell:
        """
        mv {input.fw} {output.fw}
        mv {input.rv} {output.rv}
        """

rule move_cellhash_fq:
    # Dummy rule to move the gex fastq files to the "clean" folder, in which the collapsed feature barcoding
    # fastq files will be stored.
    input:
        fw = "data/lane_merged/{sample}_CH_S1_L001_R1_001.fastq.gz",
        rv = "data/lane_merged/{sample}_CH_S1_L001_R2_001.fastq.gz"
    output:
        fw = temp("data/clean/{sample}_CH_S1_L001_R1_001.fastq.gz"),
        rv = temp("data/clean/{sample}_CH_S1_L001_R2_001.fastq.gz")
    container:
        None
    shell:
        """
        mv {input.fw} {output.fw}
        mv {input.rv} {output.rv}
        """

rule move_atac_fq:
    # Dummy rule to move the atac fastq files to the "clean" folder, in which the collapsed feature barcoding
    # fastq files will be stored.
    input:
        fw = "data/lane_merged/{sample}_ATAC_S1_L001_R1_001.fastq.gz",
        rv = "data/lane_merged/{sample}_ATAC_S1_L001_R2_001.fastq.gz",
        r3 = "data/lane_merged/{sample}_ATAC_S1_L001_R3_001.fastq.gz"
    output:
        fw = temp("data/clean/{sample}_ATAC_S1_L001_R1_001.fastq.gz"),
        rv = temp("data/clean/{sample}_ATAC_S1_L001_R2_001.fastq.gz"),
        r3 = temp("data/clean/{sample}_ATAC_S1_L001_R3_001.fastq.gz")
    container:
        None
    shell:
        """
        mv {input.fw} {output.fw}
        mv {input.rv} {output.rv}
        mv {input.r3} {output.r3}
        """