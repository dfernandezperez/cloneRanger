# No conda env used as cellranger cannot be installed in one
rule counts:
    input:
        unpack(get_sample_reads),
        bin=rules.get_cellranger.output.bin,
        genome=rules.get_reference.output.dir,
    output:
        results=directory("results/counts/{sample}/outs/raw_feature_bc_matrix"),
        mtx="results/counts/{sample}/outs/raw_feature_bc_matrix/matrix.mtx.gz",
        html=report(
            "results/counts/{sample}/outs/web_summary.html",
            caption="../reports/counts.rst",
            category="Cellranger Counts",
            subcategory="{sample}",
        ),
    params:
        introns=convert_introns(),
        n_cells=config["counts"]["n_cells"],
        mem=config["counts"]["mem"],
    log:
        "results/logs/counts/{sample}.log",
    benchmark:
        "results/benchmarks/counts/{sample}.txt"
    threads: 16
    shell:
        """
        {input.bin} \
        count \
        --nosecondary \
        {params.introns} \
        --id {wildcards.sample} \
        --transcriptome {input.genome} \
        --fastqs data \
        --sample {wildcards.sample} \
        --expect-cells {params.n_cells} \
        --localcores {threads} \
        --localmem {params.mem} \
        &> {log} && \
        mv {wildcards.sample} results/counts/{wildcards.sample}
        """
