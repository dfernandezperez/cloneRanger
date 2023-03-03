# No conda env used as cellranger cannot be installed in one
rule cellranger_count:
    input:
        lambda w: expand("data/{lane.sample_id}_S1_L00{lane.lane}_{read}_001.fastq.gz", lane=units.loc[w.sample].itertuples(), read=["R1", "R2"])
    output:
        mtx  = "results/counts/{sample}/outs/filtered_feature_bc_matrix/matrix.mtx.gz",
        html = report(
            "results/counts/{sample}/outs/web_summary.html",
            caption="../reports/counts.rst",
            category="Cellranger Counts",
            subcategory="{sample}",
        ),
    params:
        introns = convert_introns(),
        n_cells = config["cellranger_count"]["n_cells"],
        genome  = config["genome_reference"]
    log:
        "results/logs/counts/{sample}.log",
    benchmark:
        "results/benchmarks/counts/{sample}.txt"
    threads:
        RESOURCES["cellranger_count"]["cpu"]
    resources:
        mem_mb = RESOURCES["cellranger_count"]["MaxMem"]
    singularity: 
        config["cellranger_sif"]
    shell:
        """
        cellranger count \
        --nosecondary \
        {params.introns} \
        --id {wildcards.sample} \
        --transcriptome {params.genome} \
        --fastqs data \
        --sample {wildcards.sample} \
        --expect-cells {params.n_cells} \
        --localcores {threads} \
        --localmem {resources.mem_mb} \
        &> {log} && \
        # a folder in results/counts/{wildcards.sample} is automatically created due to the output declared, which 
        # is a problem to move the cellranger output files. The workaround of deleting that folder fixes that.
        rm -r results/counts/{wildcards.sample} && \
        mv {wildcards.sample} results/counts/{wildcards.sample}
        """