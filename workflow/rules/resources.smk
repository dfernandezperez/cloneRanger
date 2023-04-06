# This currently expects names as I receive them from the core
rule clean_names:
    input:
        get_fastqs,
    output:
        fw = temp("data/symlink/{sample}_{feature_bc}_S1_L00{lane}_R1_001.fastq.gz"),
        rv = temp("data/symlink/{sample}_{feature_bc}_S1_L00{lane}_R2_001.fastq.gz"),
    container:
        None
    shell:
        """
        ln -s {input[0]} {output.fw}
        ln -s {input[1]} {output.rv}
        """

rule merge_lanes:
    input:
        fw = lambda w: expand("data/symlink/{sample.sample_id}_{sample.feature_bc}_S1_L00{sample.lane}_R1_001.fastq.gz", sample=units.loc[(w.sample, w.feature_bc)].itertuples()),
        rv = lambda w: expand("data/symlink/{sample.sample_id}_{sample.feature_bc}_S1_L00{sample.lane}_R2_001.fastq.gz", sample=units.loc[(w.sample, w.feature_bc)].itertuples())
    output:
        fw = "data/lane_merged/{sample}_{feature_bc}_S1_L001_R1_001.fastq.gz",
        rv = "data/lane_merged/{sample}_{feature_bc}_S1_L001_R2_001.fastq.gz"
    log:
        "results/00_logs/merge_lanes/{sample}_{feature_bc}.log"
    shell:
        """
        cat {input.fw} > {output.fw} 2> {log}
        cat {input.rv} > {output.rv} 2>> {log}
        """

rule move_gex_fq:
    # Dummy rule to move the gex fastq files to the "clean" folder, in which the collapsed feature barcoding
    # fastq files will be stored.
    input:
        fw = "data/lane_merged/{sample}_GEX_S1_L001_R1_001.fastq.gz",
        rv = "data/lane_merged/{sample}_GEX_S1_L001_R2_001.fastq.gz"
    output:
        fw = "data/clean/{sample}_GEX_S1_L001_R1_001.fastq.gz",
        rv = "data/clean/{sample}_GEX_S1_L001_R2_001.fastq.gz"