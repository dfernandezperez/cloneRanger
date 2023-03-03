# This currently expects names as I receive them from the core
rule clean_names:
    input:
        get_fastqs,
    output:
        fw = "data/{sample}_S1_L00{lane}_R1_001.fastq.gz",
        rv = "data/{sample}_S1_L00{lane}_R2_001.fastq.gz",
    log:
        "results/logs/clean_names/{lane}_{sample}.log",
    benchmark:
        "results/benchmarks/clean_names/{lane}_{sample}.txt"
    shell:
        """
        mv {input[0]} {output.fw} &> {log}
        mv {input[1]} {output.rv} &>> {log}
        """
