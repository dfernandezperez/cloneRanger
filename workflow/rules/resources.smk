# This currently expects names as I receive them from the core
rule clean_names:
    input:
        get_fastqs,
    output:
        fw = temp("data/{sample}_{feature_bc}_{lane}_R1_symlink.fastq.gz"),
        rv = temp("data/{sample}_{feature_bc}_{lane}_R2_symlink.fastq.gz"),
    container:
        None
    shell:
        """
        ln -s {input[0]} {output.fw}
        ln -s {input[1]} {output.rv}
        """

rule merge_lanes:
    input:
        fw = lambda w: expand("data/{sample.sample_id}_{sample.feature_bc}_{sample.lane}_R1_symlink.fastq.gz", sample=units.loc[(w.sample, w.feature_bc)].itertuples()),
        rv = lambda w: expand("data/{sample.sample_id}_{sample.feature_bc}_{sample.lane}_R2_symlink.fastq.gz", sample=units.loc[(w.sample, w.feature_bc)].itertuples())
    output:
        fw = "data/{sample}_{feature_bc}_S1_L001_R1_001.fastq.gz",
        rv = "data/{sample}_{feature_bc}_S1_L001_R2_001.fastq.gz"
    log:
        "results/00_logs/merge_lanes/{sample}_{feature_bc}.log"
    message:
        "Merging fastq files from {input}"
    shell:
        """
        cat {input.fw} > {output.fw} 2> {log}
        cat {input.rv} > {output.rv} 2>> {log}
        """